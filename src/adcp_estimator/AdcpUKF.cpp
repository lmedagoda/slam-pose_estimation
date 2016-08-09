#include "AdcpUKF.hpp"
#include <base/Float.hpp>
#include <base/Logging.hpp>
#include <pose_estimation/EulerConversion.hpp>
#include <assert.h>
#include <iostream>

/** WGS-84 ellipsoid constants (Nominal Gravity Model and Earth angular velocity) **/
static const double EQUATORIAL_RADIUS = 6378137.0; /** Equatorial radius in meters **/
static const double ECC = 0.0818191908426; /** First eccentricity **/
static const double GRAVITY = 9.79766542; /** Mean value of gravity value in m/s^2 according to WGS-84 **/
static const double GRAVITY_SI = 9.80665; /** Mean value of gravity value in m/s^2 according to SI standard **/
static const double GWGS0 = 9.7803267714; /** Gravity value at the equator in m/s^2 **/
static const double GWGS1 = 0.00193185138639; /** Gravity formula constant **/
static const double EARTHW =  7.292115e-05; /** Earth angular velocity in rad/s **/ //((2.0*M_PI)/86164.0);

namespace pose_estimation
{

const std::string AdcpUKF::acceleration_measurement = "acceleration";
const std::string AdcpUKF::rotation_rate_measurement = "rotation_rate";
const std::string AdcpUKF::velocity_measurement = "velocity";
const std::string AdcpUKF::adcp_measurement = "adcp";

/** Process model for the robot orientation
 */
template <typename AdcpState>
AdcpState
processModel (const AdcpState &state, const Eigen::Vector3d& acc, const Eigen::Vector3d& omega, double delta_time,
              const AdcpUKFConfig &config, const Eigen::Vector3d& earth_rotation, const Eigen::Vector3d& gravity)
{
    AdcpState new_state(state);
    Eigen::Vector3d angular_velocity = omega - new_state.bias_gyro - new_state.orientation.inverse()*earth_rotation;
    //Eigen::Vector3d angular_velocity = omega - new_state.orientation.inverse()*earth_rotation;
    new_state.orientation.boxplus(angular_velocity, delta_time);
    
    Eigen::Vector3d acceleration = new_state.orientation*(acc - new_state.bias_acc) - gravity;
    //Eigen::Vector3d acceleration = new_state.orientation*(acc) - gravity;
    new_state.velocity.boxplus(acceleration, delta_time);
    
    new_state.position.boxplus(new_state.velocity, delta_time);
    
    Eigen::Vector3d gyro_bias_delta = (-1.0/config.gyro_bias_tau) * new_state.bias_gyro;
    new_state.bias_gyro.boxplus(gyro_bias_delta, delta_time);
    
    Eigen::Vector3d acc_bias_delta = (-1.0/config.acc_bias_tau) * new_state.bias_acc;
    new_state.bias_acc.boxplus(acc_bias_delta, delta_time);
    
    return new_state;
}

template <typename AdcpState>
Eigen::Matrix<typename AdcpState::scalar, -1, 1>
velocityMeasurementModel ( const AdcpState &state, const Eigen::Quaterniond& current_orientation )
{
    return  state.orientation.inverse() * state.velocity;
    
}

template <typename AdcpState>
Eigen::Matrix<typename AdcpState::scalar, -1, 1>
adcpMeasurementModel ( const AdcpState &state, adcp_model::ADCP_measurement_model& ADCP_model,
               AdcpUKFConfig config, boost::shared_ptr<AdcpUKF::MTK_UKF> ukf )
{
	//one beam at a time, populate a large measurement vector for the return
	ADCP_model.resetOutput();
	
	for(int i=0;i<4;i++)
	{
		Eigen::Matrix3d Rotation;
		Rotation = state.orientation.matrix();
		Eigen::Vector3d euler = Rotation.eulerAngles(0, 1, 2);
		
		boost::array<double,3> euler1; 
		euler1[0] = euler[0];
		euler1[1] = euler[1];
		euler1[2] = euler[2];
		
		boost::array<double,3> position;
		boost::array<double,3> velocity;
		
		position[0] = state.position[0];
		position[1] = state.position[1];
		position[2] = state.position[2];
		
		velocity[0] = state.velocity[0];
		velocity[1] = state.velocity[1];
		velocity[2] = state.velocity[2];		
		
		ADCP_model.setBeam(config.cell_end,position,velocity,euler1,
				   config.beam_pitch,config.beam_yaw+i*M_PI/2);
		
		bool failflag = ADCP_model.calculateInterceptsandWeightings();
		
		//std::cout << "failflag: " <<failflag << std::endl;
		
		while (unsigned int new_state = ADCP_model.isNewStates()) // this stuff probably should be handled externally. Need to guarantee through geometry that states dont get reset anyways.
		{
			if (ADCP_model.checkRemoval())
			{
				unsigned int old_state = ADCP_model.removeOldest();
				
				//remove and reinitialize the old state
				
				ukf->reset_state(old_state,0.01,0);
				ukf->reset_state(old_state+1,0.01,0);
				ukf->reset_state(old_state+2,0.01,0);
				// below emulates the process of removing the old state from the filter (zeroing it)
				
// 					active[old_state] = 0;
// 					active[old_state+1] = 0;
// 					active[old_state+2] = 0;
				
				// below emulates the process of adding the state to the filter (initializing)
				
// 					active[new_state] = 1;
// 					active[new_state+1] = 1;
// 					active[new_state+2] = 1;
			}
			
			//new states need to be activated? this is mapped to the existing covariance right?
			
// 				else
// 				{
// 					// we initialize new states, since we are filling the state vector for the first time
// 					
// 					active[new_state] = 1;
// 					active[new_state+1] = 1;
// 					active[new_state+2] = 1;
// 					
// 				}
			
		}
		
		// now calculate the predicted measurements and jacobians
		// for all valid cells as communicated by cell_start and cell_end
		ADCP_model.calculatePredictedMeasurement();
		//std::cout << "m_num:" << ADCP_model.getMeasurementNumber() << std::endl;
		for (int j=0;j<ADCP_model.getMeasurementNumber()-1;j++)
		{	
			//std::cout << " " << ADCP_model.PredictedMeasurement[j];
		}
		//std::cout << std::endl;
	}		

    //in cell order, in beam order.
    Eigen::Matrix<double, Eigen::Dynamic, 1> PredictedMeasurement;

    for (int i=0;i<ADCP_MAX_CELLS;i++)
    {
        PredictedMeasurement[4*i] = ADCP_model.PredictedMeasurement[i,0] - state.bias_adcp[0];
        PredictedMeasurement[4*i+1] = ADCP_model.PredictedMeasurement[i,1] - state.bias_adcp[1];
        PredictedMeasurement[4*i+2] = ADCP_model.PredictedMeasurement[i,2] - state.bias_adcp[2];
        PredictedMeasurement[4*i+3] = ADCP_model.PredictedMeasurement[i,3] - state.bias_adcp[3];
    }

    return  PredictedMeasurement; // may need to be trimmed. biases are just for each beam, this needs to be applied as well, but we start without them.
	
}

AdcpUKF::AdcpUKF(const AbstractFilter::FilterState& initial_state, const AdcpUKFConfig& config) : config(config)
{
    setInitialState(initial_state);
    updateFilterParamter();
}

void AdcpUKF::predictionStep(const double delta)
{
    std::map<std::string, pose_estimation::Measurement>::const_iterator acceleration = latest_measurements.find(acceleration_measurement);
    std::map<std::string, pose_estimation::Measurement>::const_iterator rotation_rate = latest_measurements.find(rotation_rate_measurement);

    if(acceleration == latest_measurements.end())
    {
        LOG_ERROR_S << "No acceleration measurement available! Skipping prediction step.";
        return;
    }

    if(rotation_rate == latest_measurements.end())
    {
        LOG_ERROR_S << "No angular velocity measurement available! Skipping prediction step.";
        return;
    }

    MTK_UKF::cov process_noise = process_noise_cov * delta; // might be pow(delta,2), depending on sensor spec.

    //process_noise.block(6,6,3,3) = 2.0 * it1->second.cov;
    // need to do uncertainty matrix calculations
    ukf->predict(boost::bind(processModel<WState>, _1, acceleration->second.mu, rotation_rate->second.mu, delta, config, earth_rotation, gravity), MTK_UKF::cov(process_noise));

}

void AdcpUKF::setFilterConfiguration(const AdcpUKFConfig& config)
{
    this->config = config;
    updateFilterParamter();
}

void AdcpUKF::correctionStepUser(const pose_estimation::Measurement& measurement)
{
    if(measurement.measurement_name == acceleration_measurement)
        latest_measurements[measurement.measurement_name] = measurement;
    else if(measurement.measurement_name == rotation_rate_measurement)
        latest_measurements[measurement.measurement_name] = measurement;
    else if(measurement.measurement_name == velocity_measurement)
    {
      std::cout << "dvl measurement:" << measurement.mu << std::endl;
        // handle velocity measurement
      std::cout << "dvl model: " << ukf->mu().orientation.inverse() * ukf->mu().velocity << std::endl;
        ukf->update(measurement.mu, boost::bind(velocityMeasurementModel<State>, _1, ukf->mu().orientation),
                        boost::bind(ukfom::id< Eigen::MatrixXd >, measurement.cov),
                        allowed_distance);
    }
    else if(measurement.measurement_name == adcp_measurement)
    {
	    std::cout << "adcp measurement:" << measurement.mu << std::endl;
	    // handle velocity measurement

        //measurement will have to be pre-ordered to match the formatting of the output measurement

	    ukf->update(measurement.mu, boost::bind(adcpMeasurementModel<State>, _1, ADCP_model, config, ukf),
			boost::bind(ukfom::id< Eigen::MatrixXd >, measurement.cov),
			allowed_distance);    
    }
    else
        LOG_ERROR_S << "Measurement " << measurement.measurement_name << " is not supported by the Adcp filter.";
}

void AdcpUKF::updateFilterParamter()
{
    double gyro_bias_var = (2.0 *pow(config.gyro_bias_std,2)) / config.gyro_bias_tau;
    double acc_bias_var = (2.0 *pow(config.acc_bias_std,2)) / config.acc_bias_tau;

    process_noise_cov = MTK_UKF::cov::Zero();
    MTK::setDiagonal(process_noise_cov, &WState::orientation, config.gyro_var*2);
    MTK::setDiagonal(process_noise_cov, &WState::velocity, config.acc_var*2);
    MTK::setDiagonal(process_noise_cov, &WState::position, 0);
    MTK::setDiagonal(process_noise_cov, &WState::bias_gyro, gyro_bias_var);
    MTK::setDiagonal(process_noise_cov, &WState::bias_acc, acc_bias_var);
    MTK::setDiagonal(process_noise_cov, &WState::bias_adcp, 0.01);
    MTK::setDiagonal(process_noise_cov, &WState::water_velocity, 0.01);

    earth_rotation[0] = cos(config.latitude)*EARTHW;
    earth_rotation[1] = 0.0;
    earth_rotation[2] = sin(config.latitude)*EARTHW;

    /** Gravity affects by the altitude (aprox the value r = EQUATORIAL_RADIUS **/
    //g = g*pow(EQUATORIAL_RADIUS/(EQUATORIAL_RADIUS+altitude), 2); //assume zero altitude for now

    gravity[0] = 0.0;
    gravity[1] = 0.0;
    gravity[2] = GWGS0*((1+GWGS1*pow(sin(config.latitude),2))/sqrt(1-pow(ECC,2)*pow(sin(config.latitude),2)));
    
    Eigen::Matrix3d Rotation;
    Rotation = ukf->mu().orientation.matrix();
    Eigen::Vector3d euler = Rotation.eulerAngles(0, 1, 2);
    boost::array<double,3> euler1; 
    euler1[0] = euler[0];
    euler1[1] = euler[1];
    euler1[2] = euler[2];
    
    boost::array<double,3> position;
    boost::array<double,3> velocity;
    
    boost::array<double,200> X_EKF; // only important to have water current states in the right spots. 200 is just a buffer size, can be made variable.
    X_EKF[0] = ukf->mu().position[0];
    X_EKF[1] = ukf->mu().position[1];
    X_EKF[2] = ukf->mu().position[2];
    
    position[0] = ukf->mu().position[0];
    position[1] = ukf->mu().position[1];
    position[2] = ukf->mu().position[2];
    
    X_EKF[3] = ukf->mu().velocity[0];
    X_EKF[4] = ukf->mu().velocity[1];
    X_EKF[5] = ukf->mu().velocity[2];
    
    velocity[0] = ukf->mu().velocity[0];
    velocity[1] = ukf->mu().velocity[1];
    velocity[2] = ukf->mu().velocity[2];
        
    X_EKF[6] = euler[0];
    X_EKF[7] = euler[1];
    X_EKF[8] = euler[2];
    
    // need to allocate the water current states
    
    ADCP_model.setConfig(config.cell_start, config.depth_cell_size, config.blank_d,
			 config.vert_grid_size,config.hori_res, X_EKF,
			 config.max_states, config.num_other_states, position, 
			 velocity, euler1);

    
}

void AdcpUKF::muToUKFState(const FilterState::Mu& mu, WState& state) const
{
    assert(mu.rows() >= WState::DOF);

    base::Orientation orientation;
    Eigen::Vector3d euler = mu.block(0, 0, 3, 1);
    EulerConversion::eulerToQuad(euler, orientation);
    state.orientation = MTK::SO3<double>(orientation);
    state.velocity = mu.block(3, 0, 3, 1);
    state.position = mu.block(6, 0, 3, 1);
    state.bias_gyro = mu.block(9, 0, 3, 1);
    state.bias_acc = mu.block(12, 0, 3, 1);
    
}

void AdcpUKF::UKFStateToMu(const WState& state, FilterState::Mu& mu) const
{
    mu.resize(WState::DOF);
    mu.setZero();

    Eigen::Vector3d euler;
    EulerConversion::quadToEuler(state.orientation, euler);
    mu.block(0, 0, 3, 1) = euler;
    mu.block(3, 0, 3, 1) = state.velocity;
    mu.block(6, 0, 3, 1) = state.position;
    mu.block(9, 0, 3, 1) = state.bias_gyro;
    mu.block(12, 0, 3, 1) = state.bias_acc;        
}

}
