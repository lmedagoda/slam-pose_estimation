#ifndef _POSE_ESTIMATION_FULL_UKF_HPP
#define _POSE_ESTIMATION_FULL_UKF_HPP

#include "FullState.hpp"
#include "FullUKFConfig.hpp"
#include <pose_estimation/Measurement.hpp>
#include <pose_estimation/UKF.hpp>
#include "adcp_measurement_model.hpp"
#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <ukfom/ukf.hpp>
#include <ukfom/mtkwrap.hpp>

namespace pose_estimation
{

class FullUKF : public UKF<FullState>
{
public:
	//typedef Manifold State;
	//typedef ukfom::mtkwrap<Manifold> WState;
	typedef ukfom::ukf<WState> MTK_UKF;
    static const std::string acceleration_measurement;
    static const std::string rotation_rate_measurement;
    static const std::string velocity_measurement;
    static const std::string adcp_measurement;

public:

    FullUKF(const FilterState& initial_state, const FullUKFConfig& config);
    virtual ~FullUKF() {}
    
    virtual void predictionStep(const double delta);

    void setFilterConfiguration(const FullUKFConfig& config);
    
    
protected:
    virtual void correctionStepUser(const pose_estimation::Measurement& measurement);

    /* This needs to called after the filter configuration was updated */
    void updateFilterParamter();

    virtual void muToUKFState(const FilterState::Mu &mu, WState& state) const;
    virtual void UKFStateToMu(const WState& state, FilterState::Mu &mu) const;

protected:
    std::map<std::string, pose_estimation::Measurement> latest_measurements;
    
    FullUKFConfig config;

    Eigen::Vector3d earth_rotation;
    Eigen::Vector3d gravity;
    adcp_model::ADCP_measurement_model ADCP_model;
};

}

#endif