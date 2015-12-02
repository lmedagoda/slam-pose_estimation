#include "PoseUKFADCP.hpp"
#include <pose_estimation/mtk_ukf/ProcessModels.hpp>
#include <pose_estimation/mtk_ukf/MeasurementModels.hpp>
#include <base/Float.hpp>
#include <base/Logging.hpp>

namespace pose_estimation
{

PoseUKF::PoseUKF() : dirty(true)
{
    last_acceleration_sample.time.microseconds = 0;

    process_noise_cov = MTK_UKF::cov::Zero();

    setDiagonal(process_noise_cov, &PoseState::position, 0.01);
    setDiagonal(process_noise_cov, &PoseState::orientation, 0.001);
    setDiagonal(process_noise_cov, &PoseState::velocity, 0.00001);
    setDiagonal(process_noise_cov, &PoseState::angular_velocity, 0.00001);
    setDiagonal(process_noise_cov, &PoseStateExtraStates::extrastates, 0.0);
}

void PoseUKF::setInitialState(const FilterState& initial_state)
{
    dirty = true;
    
    if(ukf.use_count() > 0)
    {
	ukf.reset();
    }

    WPoseState state;
    Covariance rbs_cov;
    rigidBodyStateToUKFState(body_state, state.pose_with_velocity, rbs_cov);
    MTK_UKF::cov cov = MTK_UKF::cov::Zero();
    cov.block(0,0,PoseStateExtraStates::PoseWithVelocityDOF, PoseStateExtraStates::PoseWithVelocityDOF) = rbs_cov;
    //cov.block(PoseStateExtraStates::PoseWithVelocityDOF+1, PoseStateExtraStates::PoseWithVelocityDOF+1,PoseStateExtraStates::PoseWithVelocityDOF+4,PoseStateExtraStates::PoseWithVelocityDOF+4) = 0.01*MatrixXd::Identity(4,4);
    //cov.block(PoseStateExtraStates::PoseWithVelocityDOF+5, PoseStateExtraStates::PoseWithVelocityDOF+5,PoseStateExtraStates::PoseWithVelocityDOF+64,PoseStateExtraStates::PoseWithVelocityDOF+64) = 0.1*MatrixXd::Identity(60,60);    
    MTK::setDiagonal(cov, &PoseStateExtraStates::extrastates, 0.01);
    state.extrastates.extrastates = ExtraStatesState::extrastates_type::Zero();
    
    ukf.reset(new MTK_UKF(state, cov));
}

void PoseUKF::setProcessNoiseCovariance(const FilterState::Cov& noise_cov)
{
    assert(noise_cov.rows() == process_noise_cov.rows() && noise_cov.cols() == process_noise_cov.cols());
    MTK_UKF::cov cov = MTK_UKF::cov::Zero();
    cov.block(0,0,PoseStateBias::PoseWithVelocityDOF, PoseStateBias::PoseWithVelocityDOF) = noise_cov;
    MTK::setDiagonal(cov, &PoseStateExtraStates::extrastates, 0.0);
    process_noise_cov = cov;
}

void PoseUKF::predictionStep(const double delta)
{
    if(ukf.use_count() == 0)
    {
        base::samples::RigidBodyState body_state;
        body_state.initUnknown();
 	AbstractRBSFilter::setInitialState(body_state);
    }
    dirty = true;

    if(last_acceleration_sample.time.isNull())
        ukf->predict(boost::bind(processModel<WPoseState>, _1 , delta), MTK_UKF::cov(delta * process_noise_cov));
    else
    {
        MTK_UKF::cov process_noise = process_noise_cov;
        process_noise.block(6,6,3,3) = 2.0 * last_acceleration_sample.cov_acceleration;
        ukf->predict(boost::bind(processModelWithAcceleration<WPoseState>, _1, last_acceleration_sample.acceleration, delta), MTK_UKF::cov(delta * process_noise));
    }
}

//void PoseUKF::ADCPcorrectionStep() // deal with the ADCP measurement separately since it is so different, and uses a different measurement model.

//ukf->update(sub_state, boost::bind(measurementModelADCP<WPoseState>, _1, mask),
//                     boost::bind(ukfom::id< Eigen::MatrixXd >, sub_cov),
//                     ukfom::accept_any_mahalanobis_distance<MTK_UKF::scalar_type>); // probably need an innovation gate
//     }
// }

void PoseUKF::correctionStep(const Measurement& measurement)
{
    if(ukf.use_count() == 0)
    {
        base::samples::RigidBodyState body_state;
        body_state.initUnknown();
 	AbstractRBSFilter::setInitialState(body_state);
        return;
    }

    assert(measurement.mu.rows() == MEASUREMENT_SIZE);
    assert(measurement.cov.rows() == MEASUREMENT_SIZE && measurement.cov.cols() == MEASUREMENT_SIZE);
    assert(measurement.mask.rows() == MEASUREMENT_SIZE);

    // handle acc measurements
    const Measurement::StateMask &mask = measurement.mask;
    
    if(mask[BodyStateMemberAx] > 0 && mask[BodyStateMemberAy] > 0 && mask[BodyStateMemberAz] > 0)
    {
        Eigen::Vector3d acc = measurement.mu.block(12,0,3,1);
        Eigen::Matrix3d acc_cov = measurement.cov.block(12,12,3,3);
        if(base::isnotnan(acc_cov) && base::isnotnan(acc))
        {
            last_acceleration_sample.acceleration = acc;
            last_acceleration_sample.cov_acceleration = acc_cov;
            last_acceleration_sample.time = measurement.time;
        }
        else
            LOG_ERROR("Acceleration covariance or acceleration sample contains NaN values!");
    }

    // handle body state measurements
    unsigned body_state_members = 0;
    for(unsigned i = 0; i < WPoseState::DOF; i++)
        body_state_members += mask[i];
    if(body_state_members > 0)
    {
        dirty = true;

        WPoseState state;
        muToUKFState(measurement.mu, state);
        MTK_UKF::cov cov = measurement.cov.block(0,0,WPoseState::DOF,WPoseState::DOF);

        // check state for NaN values
        Eigen::Matrix<unsigned, WPoseState::DOF, 1> sub_mask = mask.block(0,0,WPoseState::DOF,1);
        Eigen::Matrix<WPoseState::scalar_type, WPoseState::DOF, 1> state_vector = state.getStateVector();
        for(unsigned i = 0; i < WPoseState::DOF; i++)
        {
            if(sub_mask(i,0) != 0 && base::isNaN<WPoseState::scalar_type>(state_vector(i,0)))
            {
                // handle NaN values in state
                LOG_ERROR("State contains NaN values! This Measurement will be excluded from correction step.");
                sub_mask(i,0) = 0;
            }
        }

        // check covariance matrix for NaN values
        for(unsigned i = 0; i < WPoseState::DOF; i++)
            for(unsigned j = 0; j < WPoseState::DOF; j++)
            {
                if(sub_mask(i,0) != 0 && sub_mask(j,0) != 0 && base::isNaN<Covariance::Scalar>(cov(i,j)))
                {
                    // handle NaN variances
                    LOG_ERROR("Covariance contains NaN values! This Measurement will be skipped.");
                    return;
                }
                else if(i==j && cov(i,j) == 0.0)
                {
                    // handle zero variances
                    LOG_WARN("Covariance diagonal contains zero values. Override them with %d", 1e-9);
                    cov(i,j) = 1e-9;
                }
            }

        Eigen::Matrix<WPoseState::scalar, -1, 1> sub_state = state.getSubStateVector(sub_mask);
        unsigned measurement_size = sub_state.rows();

        std::vector<unsigned> m_index2mask_index;
        for(unsigned i = 0; i < WPoseState::DOF; i++)
        {
            if(sub_mask(i) > 0)
            {
                m_index2mask_index.push_back(i);
            }
        }
        Eigen::Matrix<WPoseState::scalar, -1, -1> sub_cov(measurement_size, measurement_size);
        sub_cov.setZero();
        for(unsigned i = 0; i < measurement_size; i++)
        {
            for(unsigned j = 0; j < measurement_size; j++)
            {
                sub_cov(i,j) = cov(m_index2mask_index[i], m_index2mask_index[j]);
            }
        }

        // apply new measurement
        ukf->update(sub_state, boost::bind(measurementModel<WPoseState>, _1, sub_mask),
                    boost::bind(ukfom::id< Eigen::MatrixXd >, sub_cov),
                    ukfom::accept_any_mahalanobis_distance<MTK_UKF::scalar_type>);
    }
}

const AbstractFilter::FilterState& PoseUKF::getCurrentState()
{
    if(dirty)
    {
        UKFStateToMu(ukf->mu(), state.mu);
        state.cov = ukf->sigma();

        dirty = false;
    }
    
//+    {
//+        PoseState::cov cov = ukf->sigma().block(0, 0, PoseState::DOF, PoseState::DOF);
//+	UKFStateToRigidBodyState(ukf->mu().pose_with_velocity, cov, body_state);
//+    }
    
    return state;
}

//+base::VectorXd PoseUKF::getFullState()
//+{
//+    PoseUKF::WPoseState mu = ukf->mu();
//+    Eigen::Matrix<WPoseState::scalar_type, PoseState::DOF, 1> pose_state = mu.getStateVector();
//+    base::VectorXd state(WPoseState::DOF);
//+    state.block(WPoseState::PoseWithVelocityIdx,0,WPoseState::PoseWithVelocityDOF,1) = pose_state;
//+    state.block(WPoseState::BiasIdx,0,WPoseState::BiasDOF,1) = mu.bias.bias;
//+    return state;
//+}
//+
//+base::MatrixXd PoseUKF::getFullCovariance()
//+{
//+    MTK_UKF::cov sigma = ukf->sigma();
//+    base::MatrixXd cov(sigma.rows(), sigma.rows());
//+    cov = sigma;
//+    return cov;
//+}

void PoseUKF::muToUKFState(const StateAndCovariance::Mu& mu, PoseUKF::WPoseState& state)
{
    assert(mu.rows() >= PoseUKF::WPoseState::DOF);

    state.position = mu.block(0, 0, 3, 1);
    base::Orientation orientation;
    Eigen::Vector3d euler = mu.block(3, 0, 3, 1);
    eulerToQuad(euler, orientation);
    state.orientation = MTK::SO3<double>(orientation);
    state.velocity = mu.block(6, 0, 3, 1);
    Eigen::Vector3d angle_axis;
    Eigen::Vector3d euler_velocity = mu.block(9, 0, 3, 1);
    eulerAngleVelocityToAngleAxis(euler_velocity, angle_axis);
    state.angular_velocity = angle_axis;

}

void PoseUKF::UKFStateToMu(const PoseUKF::WPoseState& state, StateAndCovariance::Mu& mu)
{
    mu.resize(PoseUKF::WPoseState::DOF);
    mu.setZero();

    mu.block(0, 0, 3, 1) = state.position;
    Eigen::Vector3d euler;
    base::Orientation orientation = state.orientation;
    quadToEuler(orientation, euler);
    mu.block(3, 0, 3, 1) = euler;
    mu.block(6, 0, 3, 1) = state.velocity;
    Eigen::Vector3d euler_velocity;
    Eigen::Vector3d angular_velocity = state.angular_velocity;
    angleAxisToEulerAngleVelocity(angular_velocity, euler_velocity);
    mu.block(9, 0, 3, 1) = euler_velocity;
}

//void PoseUKF::rigidBodyStateToUKFState(const base::samples::RigidBodyState& body_state, PoseUKF::PoseState& state, PoseState::cov& covariance)
// {
//     state.position = body_state.position;
//    state.orientation = MTK::SO3<double>(body_state.orientation);
//     state.velocity = body_state.velocity;
//     state.angular_velocity = body_state.angular_velocity;
//
//     covariance.setZero();
//     covariance.block(0, 0, 3, 3) = body_state.cov_position;
//     covariance.block(3, 3, 3, 3) = body_state.cov_orientation;
//     covariance.block(6, 6, 3, 3) = body_state.cov_velocity;
//     covariance.block(9, 9, 3, 3) = body_state.cov_angular_velocity;
// }

//void PoseUKF::UKFStateToRigidBodyState(const PoseUKF::PoseState& state, const PoseState::cov& covariance, base::samples::RigidBodyState& body_state)
// {
//     body_state.position = state.position;
//     body_state.orientation = state.orientation;
//     body_state.velocity = body_state.orientation * state.velocity;
//     body_state.angular_velocity = state.angular_velocity;
//     
//     body_state.cov_position = covariance.block(0, 0, 3, 3);
//     body_state.cov_orientation = covariance.block(3, 3, 3, 3);
//     body_state.cov_velocity = covariance.block(6, 6, 3, 3);
//     body_state.cov_angular_velocity = covariance.block(9, 9, 3, 3);
// }
}
