#ifndef _POSE_ESTIMATION_FULL_UKF_CONFIG_HPP
#define _POSE_ESTIMATION_FULL_UKF_CONFIG_HPP

namespace pose_estimation
{

struct FullUKFConfig
{
    double gyro_var; // (rad/s)/sqrt(Hz)
    double acc_var; // (m/s)/sqrt(Hz)
    double gyro_bias_std; // rad/s
    double acc_bias_std;  // m/s^2
    double gyro_bias_tau; // in seconds
    double acc_bias_tau;
    double latitude; // in radians

    double cell_start;//default = 1;
    double cell_end;//default = 10;
    double depth_cell_size;//default = 1;
    double blank_d;//default = 0.1;
    double vert_grid_size;//default = 15;
    double hori_res;//default = 50;
    double beam_pitch;//default = M_PI/6;
    double beam_yaw;//default = M_PI/4;
    bool failflag;//default = false;
    double max_states;//default = 76;
    double num_other_states;//default = 16;
    unsigned short new_state;
    unsigned short old_state;
    

	
};



}

#endif