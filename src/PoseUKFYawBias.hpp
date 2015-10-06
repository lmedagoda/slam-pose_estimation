#ifndef _POSE_ESTIMATION_POSE_UKF_YB_HPP
#define _POSE_ESTIMATION_POSE_UKF_YB_HPP

#include <pose_estimation/Measurement.hpp>
#include <pose_estimation/AbstractFilter.hpp>
#include <pose_estimation/mtk_ukf/PoseWithVelocity.hpp>
#include <iostream>
#include <ukfom/ukf.hpp>
#include <ukfom/mtkwrap.hpp>
#include <boost/shared_ptr.hpp>

namespace pose_estimation
{
    
class PoseUKFYawBias : public AbstractFilter
{
public:
    PoseUKFYawBias();
    
    virtual void setInitialState(const base::samples::RigidBodyState &body_state);
    virtual void setProcessNoiseCovariance(const Covariance& noise_cov);
    virtual void predictionStep(const double delta);
    virtual void correctionStep(const Measurement &measurement);
    virtual const base::samples::RigidBodyState& getCurrentState();
    
protected:
    typedef MTK::SO3<double> MTKRotationType;
    typedef PoseWithVelocity<MTKRotationType> PoseState;
    typedef MTK::vect<1, double> YawBias;
    //Define a Compound manifold
    MTK_BUILD_MANIFOLD(PoseStateADCP, 
        ((PoseState, pose))      
        ((YawBias, yaw_bias))    
    );                     
    typedef ukfom::mtkwrap<PoseStateADCP> WPoseState;
    typedef ukfom::ukf<WPoseState> MTK_UKF;
    
    void rigidBodyStateToUKFState(const base::samples::RigidBodyState &body_state, WPoseState& state, MTK_UKF::cov& covariance);
    void UKFStateToRigidBodyState(const WPoseState& state, const MTK_UKF::cov& covariance, base::samples::RigidBodyState &body_state);
    
protected:    
    boost::shared_ptr<MTK_UKF> ukf;
    bool dirty;
    base::samples::RigidBodyState body_state;
    base::samples::RigidBodyAcceleration last_acceleration_sample;
    Covariance process_noise_cov;
};

}

#endif
