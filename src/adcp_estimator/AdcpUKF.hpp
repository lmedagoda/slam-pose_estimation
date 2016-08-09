#ifndef _POSE_ESTIMATION_ADCP_UKF_HPP
#define _POSE_ESTIMATION_ADCP_UKF_HPP

#include "AdcpState.hpp"
#include "AdcpUKFConfig.hpp"
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

class AdcpUKF : public UKF<AdcpState>
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

    AdcpUKF(const FilterState& initial_state, const AdcpUKFConfig& config);
    virtual ~AdcpUKF() {}
    
    virtual void predictionStep(const double delta);

    void setFilterConfiguration(const AdcpUKFConfig& config);
    
    
protected:
    virtual void correctionStepUser(const pose_estimation::Measurement& measurement);

    /* This needs to called after the filter configuration was updated */
    void updateFilterParamter();

    virtual void muToUKFState(const FilterState::Mu &mu, WState& state) const;
    virtual void UKFStateToMu(const WState& state, FilterState::Mu &mu) const;

protected:
    std::map<std::string, pose_estimation::Measurement> latest_measurements;
    
    AdcpUKFConfig config;

    Eigen::Vector3d earth_rotation;
    Eigen::Vector3d gravity;
    adcp_model::ADCP_measurement_model ADCP_model;
};

}

#endif
