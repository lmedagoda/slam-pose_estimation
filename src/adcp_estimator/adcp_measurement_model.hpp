#ifndef ADCP_MEASUREMENT_MODEL_H
#define ADCP_MEASUREMENT_MODEL_H

// these definitions cannot be changed without remaking the matlab code. Present iteration is not flexible to changes due to array indexing from MATLAB.
#define ADCP_MAX_INTERCEPTS 7
#define ADCP_MAX_CELLS 30
#define ADCP_VERTEX_NUMBER 8
#define ADCP_MAX_STATES 200

#include <iostream>
#include <math.h>
#include "rtwtypes.h"
#include <boost/array.hpp>
#include <vector>
#include <map>
#include <string.h>
#include <boost/circular_buffer.hpp>
#include <eigen3/Eigen/Sparse>
namespace adcp_model
{
	class ADCP_measurement_model
	{
	public:
	ADCP_measurement_model();


	~ADCP_measurement_model();
	void setBeam(double cell_end_,
			boost::array<double,3>& position_,
			boost::array<double,3>& velocity_,
			boost::array<double,3>& euler_,
			double beam_pitch_,
			double beam_yaw_);

	void resetOutput(void);
	
	void setConfig(const unsigned short cell_start_,
		       const double depth_cell_size_, // size of the ADCP measurement cells
		const double blank_d_,
		const double vert_grid_size_, // size of the depth cells in the state vector (assumed to be equal or larger than DepthCellSize);
	const double hori_res_,
	const boost::array<double,200>& X_EKF_,
	const unsigned short max_states_,
	const unsigned short num_other_states_,
	const boost::array<double,3>& position_,
	const boost::array<double,3>& velocity_,
	const boost::array<double,3>& euler_); // initialize constants that don't change from beam to beam, and measurement to measurement

	// 1)
	void initializeModel(void); // All the arrays used are zeroed, since we want to use them for checking later.
	bool calculateInterceptsandWeightings(void); // returns a failflag, which is true if there was a failure.

	// 2) [external signal] Compare observed states with present states. Initialize new ones if necessary. If present state number exceeds maximum, marginalize out the oldest state.
	// will call setLocalMap() since this can be done in parallel

	unsigned short isNewStates(void); //returns non-zero if there are new states to add to the filter due to new observing of cells. States are added according to a ring.
	//The return value is the location of the new state.
	bool checkRemoval(void);
	unsigned short removeOldest(void); // also returns the location of the oldest state, so that it can be replaced.

	// 3) [external signal] Link new states to existing states with Kalman updates (neighbourhood correlation). This is optional for now and is not initially implemented.

	// 4) h_x_ADCP_dfki
	bool calculatePredictedMeasurement(void); // returns the fail flag. If true then there was a problem.
	unsigned short getMeasurementNumber(void);

	// Also calculates Jacobians

	Eigen::Matrix<double,ADCP_MAX_CELLS*4,1> PredictedMeasurement; // *4 because this can store the all the beams.
	//    Eigen::Matrix<double,ADCP_MAX_CELLS*4,ADCP_MAX_STATES> Jacobian; // you can easily convert this to other types through referencing
	Eigen::SparseMatrix<double> Jacobian; // use triplets to construct this.
	// you can easily convert this to other types through referencing
	// Eigen supports using this through Map: http://eigen.tuxfamily.org/dox-devel/group__TutorialMapClass.html
	// probably just natively use Eigen?

	void test(void);

	private:

	unsigned short measurement_counter;

	double current_weightings[ADCP_VERTEX_NUMBER*ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS]; // 8 weights, 7 possible intercepts, 30 measurement cells maximum

	int grid_loc_row[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS], grid_loc_col[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS], grid_loc_E[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS];
	unsigned int active_states[ADCP_MAX_CELLS]; // tracks the number of intercepts for each cell

	boost::array<double,200> X_EKF; // make sure you get the referencing right. All references are from the first ADCP state
	// This should be only the portion of the state vector which are water current velocities

	unsigned short oldest_location, newest_location; // since the filter treats the state vector as a ring buffer for water current states, we track the oldest and newest states in that ring buffer
	// oldest location increments when a state is removed. newest location increments when a state is added. Remember, this increments by 3 since each state has XYZ velocity.
	// initialize these as 1?
	bool state_filled;

	std::map<boost::array<int,3>,unsigned short> Map; // key is xyz coordinate. index is the state. We know the oldest and youngest states by tracking them. This will get passed to h_x_ADCP_dfki.
	// allows ~constant time element access of keys

	//    boost::circular_buffer<boost::array<int,3>> xyzList(ADCP_MAX_STATES);
	//    boost::circular_buffer<(boost::array<int,3>)> xyzList;
	std::deque<boost::array<int,3> > xyzList; // x,y,z, in order of intialization. Will be used when removal and addition of states is invoked, to match keys to index in map.
	// allows constant time removal and addition of states. Used to find the map key for removal as well.

	//    static const double P_init[3] = {0,0,0}; // this is to preserve the relationship with the MATLAB code, which does not have the same indexing on the map. Should do it within the function.

	real_T cell_start;
	real_T cell_end;
	boost::array<double,3> position;
	boost::array<double,3> velocity;
	boost::array<double,3> euler;
	real_T depth_cell_size; // size of the ADCP measurement cells
	real_T blank_d;
	real_T vert_grid_size; // size of the depth cells in the state vector (assumed to be equal or larger than DepthCellSize);
	real_T hori_res;
	real_T beam_pitch;
	real_T beam_yaw;
	bool failflag;
	real_T  max_states; // maximum amount of water current states that are tracked (initially set this to 60)
	real_T  num_other_states;

	double J_vc[200]; // for the jacobian calculation

	double dz_dh[3], dz_dv[3];

	//    setLocalMap(); // takes in the global hashtable map, and creates local map (3D array with states).. is this too slow? Happens during the check for observed states and present states
	// may not be necessary, we can just reference the Map directly.

	// for tracking the state checking (based on what we are up to in checking grid_loc_row etc.
	// used by isNewStates
	unsigned short new_state_track, new_state_track_cell;

	//1) [current_weightings,grid_loc_row,grid_loc_col,grid_loc_E] = findInterceptsandWeightings(position,...
	//     euler,depth_cell_size,cell_start,cell_end,blank_d,hori_res,beam_pitch,beam_yaw,P_init,vert_grid_size)
	void findInterceptsandWeightings(void);

	//    void findInterceptsandWeightings(const real_T position[3], const real_T euler[3], real_T depth_cell_size, real_T cell_start,
	//    real_T cell_end, real_T blank_d, real_T hori_res, real_T beam_pitch, real_T beam_yaw, const real_T P_init[3], real_T vert_grid_size,
	//    real_T current_weightings[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS*ADCP_VERTEX_NUMBER], real_T grid_loc_row[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS],
	//    real_T grid_loc_col[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS], real_T grid_loc_E[ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS], bool *failflag);


	// this happens for each beam

	//1) a) Intercepts [i,i_num,failflag] = find_intercepts_3d_4dof_dfki(x1,y1,z1,x2,y2,z2,hori_res,depth_cell_size)
	// one beam at a time
	void find_intercepts_3d_4dof_dfki(real_T x1, real_T b_y1, real_T z1, real_T x2,
					real_T y2, real_T z2, real_T hori_res, real_T depth_cell_size, real_T i[21],
	real_T *i_num, bool *failflag); // be careful when there are single values since they don't get treated as pointers in the function (like failflag)

	void find_weighting_trilin(real_T x0, real_T b_y0, real_T z0, real_T xe, real_T ye, real_T ze, real_T l1, real_T l2, real_T w[8]);

	void eml_sort(const real_T x_data[5], const int32_T x_size[1], real_T
	y_data[5], int32_T y_size[1], int32_T idx_data[5], int32_T
	idx_size[1]);

	// 4) h_x_ADCP requires the local map
	double h_x_ADCP_dfki(uint16_T cell_index);

	// 5) Jacobians with respect to water current velocities (one measurement at a time)
	void calculateWaterCurrentVelocityJacobian(uint16_T cell_index);

	};
}

#endif // ADCP_MEASUREMENT_MODEL_H
