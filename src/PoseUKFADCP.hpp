#ifndef _POSE_ESTIMATION_POSE_UKF_ADCP_HPP
#define _POSE_ESTIMATION_POSE_UKF_ADCP_HPP

#include <pose_estimation/Measurement.hpp>
#include <pose_estimation/AbstractFilter.hpp>
#include <pose_estimation/mtk_ukf/ExtraStates.hpp>
#include <pose_estimation/mtk_ukf/PoseWithVelocityExtraStates.hpp>
#include <pose_estimation/mtk_ukf/PoseWithVelocity.hpp>
#include <iostream>
#include <ukfom/ukf.hpp>
#include <ukfom/mtkwrap.hpp>
#include <boost/shared_ptr.hpp>

namespace pose_estimation
{
    
class PoseUKF : public AbstractRBSFilter
{
public:
    PoseUKF();
    virtual ~PoseUKF() {}
    
    virtual void setInitialState(const FilterState& initial_state);
    virtual void setProcessNoiseCovariance(const FilterState::Cov& noise_cov);
    virtual void predictionStep(const double delta);
    virtual void correctionStep(const Measurement& measurement);
    virtual const FilterState& getCurrentState();
    //virtual base::VectorXd getFullState();
    //virtual base::MatrixXd getFullCovariance();
    
protected:
    typedef MTK::SO3<double> MTKRotationType;
    typedef PoseWithVelocity<MTKRotationType> PoseState;
    typedef ExtraStates<64> ExtraStatesState; // 4 beam bias states + 20 cells with 3 dimensions each = 64
    typedef PoseWithVelocityExtraStates<PoseState, ExtraStatesState> PoseStateExtraStates;
    typedef ukfom::mtkwrap<PoseStateExtraStates> WPoseState;
    typedef ukfom::ukf<WPoseState> MTK_UKF;

    static void muToUKFState(const FilterState::Mu &mu, WPoseState& state);
    static void UKFStateToMu(const WPoseState& state, FilterState::Mu &mu);
    
//+    void rigidBodyStateToUKFState(const base::samples::RigidBodyState &body_state, PoseState& state, PoseState::cov& covariance);
//+    void UKFStateToRigidBodyState(const PoseState& state, const PoseState::cov& covariance, base::samples::RigidBodyState &body_state);
//+
//+    template<class Base, class T, int idx, int cov_size>
//+    void setDiagonal(Eigen::Matrix<typename Base::scalar, cov_size, cov_size> &cov,
//+                    MTK::SubManifold<T, idx> Base::*, const typename Base::scalar &val)
//+    {
//+        cov.diagonal().template segment<T::DOF>(idx).setConstant(val);
//+    }
    

protected:
    boost::shared_ptr<MTK_UKF> ukf;
    bool dirty;
    FilterState state;
    base::samples::RigidBodyAcceleration last_acceleration_sample;
    MTK_UKF::cov process_noise_cov;
};



}

#endif
