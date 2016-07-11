#include "adcp_measurement_model.hpp"

using namespace adcp_model;

ADCP_measurement_model::~ADCP_measurement_model()
{

}

ADCP_measurement_model::ADCP_measurement_model()
{
    // zero all the arrays which will be used. Make a function to do this again. Or will that function always do the zeroing in the first calculate step?
    // current_weightings
    // grid_loc_row

    oldest_location = num_other_states+1;
    newest_location = num_other_states+1;
    state_filled = false;

}

void ADCP_measurement_model::setConfig(const unsigned short cell_start_,
	       const double depth_cell_size_, // size of the ADCP measurement cells
	       const double blank_d_,
	       const double vert_grid_size_, // size of the depth cells in the state vector (assumed to be equal or larger than DepthCellSize);
const double hori_res_,
const boost::array<double,200>& X_EKF_,
const unsigned short max_states_,
const unsigned short num_other_states_,
const boost::array<double,3>& position_,
const boost::array<double,3>& velocity_,
const boost::array<double,3>& euler_)
{
	cell_start = cell_start_;
	depth_cell_size = depth_cell_size_; 
	blank_d = blank_d_;
	vert_grid_size = vert_grid_size_;
	hori_res = hori_res_;
	X_EKF = X_EKF_;
	max_states = max_states_;
	num_other_states = num_other_states_;
	position = position_;
	velocity = velocity_;
	euler = euler_;
}

// const unsigned short cell_start_,
// const double depth_cell_size_, // size of the ADCP measurement cells
// const double blank_d_,
// const double vert_grid_size_, // size of the depth cells in the state vector (assumed to be equal or larger than DepthCellSize);
// const double hori_res_,
// const boost::array<double,200>& X_EKF_,
// const unsigned short max_states_,
// const unsigned short num_other_states_,
// boost::array<double,3>& position_,
// boost::array<double,3>& velocity_,
// boost::array<double,3>& euler_) : // initialize constants that don't change from beam to beam, and measurement to measurement


void ADCP_measurement_model::resetOutput(void)
{
    unsigned int i;//,j;

    for (i=0;i<ADCP_MAX_CELLS*4;i++)
    {
        PredictedMeasurement[i]=0;
//        for (j=0;j<ADCP_MAX_STATES;j++)
//        {
//            Jacobian(i,j)=0;
//        }
    }
    measurement_counter = 0;

}

void ADCP_measurement_model::setBeam(double cell_end_,
                                     boost::array<double,3>& position_,
                                     boost::array<double,3>& velocity_,
                                     boost::array<double,3>& euler_,
                                     double beam_pitch_,
                                     double beam_yaw_)
{
    cell_end = cell_end_;
    position = position_;
    velocity = velocity_;
    euler = euler_;
    beam_pitch = beam_pitch_;
    beam_yaw = beam_yaw_;
}

void ADCP_measurement_model::initializeModel(void)
{
    unsigned short i;
    // zero everything that needs to be re-tracked
    new_state_track_cell = 0;

    for(i=0;i<ADCP_MAX_CELLS;i++)
    {
        active_states[i]=0;
    }
}

bool ADCP_measurement_model::calculateInterceptsandWeightings(void)
{
//    unsigned int i;

    initializeModel();

    findInterceptsandWeightings();

    //    for (i=0;i<ADCP_VERTEX_NUMBER*ADCP_MAX_INTERCEPTS*ADCP_MAX_CELLS;i++)
    //    {
    //        std::cout << current_weightings[i];
    //        std::cout << " ";
    //    }

    return failflag;
}

unsigned short ADCP_measurement_model::isNewStates(void)
{
    boost::array<int,3> ref;
    unsigned short i = 0;
    static const int8_T gmt[8][3] = {{0,0,0},
                                     {1,0,0},
                                     {0,1,0},
                                     {0,0,1},
                                     {1,0,1},
                                     {0,1,1},
                                     {1,1,0},
                                     {1,1,1}};

    // compare the existing states in the map to the observed states
    while (new_state_track_cell <= cell_end)
    {
        new_state_track = 0;
        while (new_state_track < ADCP_MAX_INTERCEPTS)
        {
            if (new_state_track<active_states[new_state_track_cell]) // we have a new state perhaps
            {
                for (i=0;i<8;i++)
                {
                    ref[0] = grid_loc_row[new_state_track + ADCP_MAX_INTERCEPTS*new_state_track_cell] + gmt[i][1];
                    ref[1] = grid_loc_col[new_state_track + ADCP_MAX_INTERCEPTS*new_state_track_cell] + gmt[i][2];
                    ref[2] = grid_loc_E[new_state_track + ADCP_MAX_INTERCEPTS*new_state_track_cell] + gmt[i][3];

                    if (Map[ref]==0) // check if its in the map
                    {

                        if (newest_location+3 > max_states)
                        {
                            newest_location = num_other_states+1;
                            if (state_filled==false)
                            {
                                state_filled = true;
                            }
                        }

                        Map[ref] = newest_location;
                        xyzList.push_back(ref);
                        std::cout << "Key: ";
                        std::cout << ref[0] << " " << ref[1] << " " << ref[2] << std::endl;
                        std::cout << "Newest added: ";
                        std::cout << Map[ref] << std::endl;
                        newest_location += 3;
                        return newest_location-3;
                    }
                }

            }


            new_state_track++;

        }
        new_state_track_cell++;
    }

    return 0;
}

bool ADCP_measurement_model::checkRemoval(void)
{
    if ((newest_location-3 == oldest_location)&&(state_filled))
    {
        return true;
    }
    else
    {
        return false;
    }
}

unsigned short ADCP_measurement_model::removeOldest(void)
{
    // find the oldest state, remove it from the map and list, and relay the position of the oldest state to the filter so that it can be removed there too
    boost::array<int,3> ref;
    unsigned short location;

    ref = xyzList.front();
    xyzList.pop_front();
    location = Map[ref];

    std::cout << "Oldest removed: ";
    std::cout << Map[ref] << std::endl;

    Map.erase(ref);

    oldest_location += 3;
    if (oldest_location+3>max_states)
    {
        oldest_location = num_other_states+1;
    }
    return location;
}

bool ADCP_measurement_model::calculatePredictedMeasurement(void)
{
    unsigned short cell;
    bool fail = false;
    // calculate the predicted measurement, and the accompanying Jacobian, and place them in the output arrays

    // we can do this for the present beam
    for(cell=cell_start;cell<cell_end;cell++)
    {
        PredictedMeasurement[measurement_counter] = h_x_ADCP_dfki(cell);
        if ((fail == false)&&(failflag==true))
        {
            fail = true;
        }
        measurement_counter++;
    }

    return fail;
}

unsigned short ADCP_measurement_model::getMeasurementNumber(void)
{
    return measurement_counter;
}

/*
               * function [i,i_num,failflag] = find_intercepts_3d_4dof_dfki(x1,y1,z1,x2,y2,z2,hori_res,depth_cell_size)
               */
void ADCP_measurement_model::find_intercepts_3d_4dof_dfki(real_T x1, real_T b_y1, real_T z1, real_T x2, real_T y2, real_T z2, real_T hori_res,
                                                          real_T depth_cell_size, real_T i[21], real_T *i_num, bool *failflag)
{
    real_T y;

    /*  edit history: */
    /*  5/4/11 - created, refer to notes for references frames */
    /*  x1 < x2 so we can always work left to right */
    /* 'find_intercepts_3d_4dof_dfki:6' i = zeros(7,3); */
    memset(&i[0], 0, 21U * sizeof(real_T));

    /* 'find_intercepts_3d_4dof_dfki:7' i_num = 1; */
    *i_num = 1.0;

    /* 'find_intercepts_3d_4dof_dfki:8' failflag = 0; */
    *failflag = 0.0;

    /* 'find_intercepts_3d_4dof_dfki:10' if x1 > x2 */
    if (x1 > x2) {
        /* 'find_intercepts_3d_4dof_dfki:11' failflag = 1; */
        *failflag = 1.0;
    }

    /*  else */
    /*  depth_cell intercepts first */
    /* 'find_intercepts_3d_4dof_dfki:17' if y2 > y1 */
    if (y2 > b_y1) {
        /*  slant down to the right */
        /* 'find_intercepts_3d_4dof_dfki:19' if floor(y2/depth_cell_size) == ceil(y1/depth_cell_size) */
        if (floor(y2 / depth_cell_size) == ceil(b_y1 / depth_cell_size)) {
            /*  1 depth_cell intercept */
            /* 'find_intercepts_3d_4dof_dfki:21' y = floor(y2/depth_cell_size)*depth_cell_size; */
            y = floor(y2 / depth_cell_size) * depth_cell_size;

            /* 'find_intercepts_3d_4dof_dfki:22' x = (x1-x2)*(y-y1)/(y1-y2) + x1; */
            /* 'find_intercepts_3d_4dof_dfki:23' z = (z1-z2)*(y-y1)/(y1-y2) + z1; */
            /* 'find_intercepts_3d_4dof_dfki:24' i(i_num,1) = x; */
            i[0] = (x1 - x2) * (y - b_y1) / (b_y1 - y2) + x1;

            /* 'find_intercepts_3d_4dof_dfki:25' i(i_num,2) = y; */
            i[7] = y;

            /* 'find_intercepts_3d_4dof_dfki:26' i(i_num,3) = z; */
            i[14] = (z1 - z2) * (y - b_y1) / (b_y1 - y2) + z1;

            /* 'find_intercepts_3d_4dof_dfki:27' i_num = i_num + 1; */
            *i_num = 2.0;
        } else if (floor(y2 / depth_cell_size) == ceil(b_y1 / depth_cell_size) + 1.0)
        {
            /* 'find_intercepts_3d_4dof_dfki:28' elseif floor(y2/depth_cell_size) == ceil(y1/depth_cell_size)+1 */
            /*  2 depth_cell intercepts */
            /* 'find_intercepts_3d_4dof_dfki:30' y = floor(y2/depth_cell_size)*depth_cell_size; */
            y = floor(y2 / depth_cell_size) * depth_cell_size;

            /* 'find_intercepts_3d_4dof_dfki:31' x = (x1-x2)*(y-y1)/(y1-y2) + x1; */
            /* 'find_intercepts_3d_4dof_dfki:32' z = (z1-z2)*(y-y1)/(y1-y2) + z1; */
            /* 'find_intercepts_3d_4dof_dfki:33' i(i_num,1) = x; */
            i[0] = (x1 - x2) * (y - b_y1) / (b_y1 - y2) + x1;

            /* 'find_intercepts_3d_4dof_dfki:34' i(i_num,2) = y; */
            i[7] = y;

            /* 'find_intercepts_3d_4dof_dfki:35' i(i_num,3) = z; */
            i[14] = (z1 - z2) * (y - b_y1) / (b_y1 - y2) + z1;

            /* 'find_intercepts_3d_4dof_dfki:36' i_num = i_num + 1; */
            /* 'find_intercepts_3d_4dof_dfki:38' y = ceil(y1/depth_cell_size)*depth_cell_size; */
            y = ceil(b_y1 / depth_cell_size) * depth_cell_size;

            /* 'find_intercepts_3d_4dof_dfki:39' x = (x1-x2)*(y-y1)/(y1-y2) + x1; */
            /* 'find_intercepts_3d_4dof_dfki:40' z = (z1-z2)*(y-y1)/(y1-y2) + z1; */
            /* 'find_intercepts_3d_4dof_dfki:41' i(i_num,1) = x; */
            i[1] = (x1 - x2) * (y - b_y1) / (b_y1 - y2) + x1;

            /* 'find_intercepts_3d_4dof_dfki:42' i(i_num,2) = y; */
            i[8] = y;

            /* 'find_intercepts_3d_4dof_dfki:43' i(i_num,3) = z; */
            i[15] = (z1 - z2) * (y - b_y1) / (b_y1 - y2) + z1;

            /* 'find_intercepts_3d_4dof_dfki:44' i_num = i_num + 1; */
            *i_num = 3.0;
        } else {
            /* 'find_intercepts_3d_4dof_dfki:45' else */
            /*  no depth_cell_intercepts */
        }
    } else {
        /* 'find_intercepts_3d_4dof_dfki:48' else */
        /*  slant up to the right */
        /* 'find_intercepts_3d_4dof_dfki:49' if floor(y1/depth_cell_size) == ceil(y2/depth_cell_size) */
        if (floor(b_y1 / depth_cell_size) == ceil(y2 / depth_cell_size)) {
            /*  1 depth_cell intercept */
            /* 'find_intercepts_3d_4dof_dfki:51' y = floor(y1/depth_cell_size)*depth_cell_size; */
            y = floor(b_y1 / depth_cell_size) * depth_cell_size;

            /* 'find_intercepts_3d_4dof_dfki:52' x = (x1-x2)*(y-y1)/(y1-y2) + x1; */
            /* 'find_intercepts_3d_4dof_dfki:53' z = (z1-z2)*(y-y1)/(y1-y2) + z1; */
            /* 'find_intercepts_3d_4dof_dfki:54' i(i_num,1) = x; */
            i[0] = (x1 - x2) * (y - b_y1) / (b_y1 - y2) + x1;

            /* 'find_intercepts_3d_4dof_dfki:55' i(i_num,2) = y; */
            i[7] = y;

            /* 'find_intercepts_3d_4dof_dfki:56' i(i_num,3) = z; */
            i[14] = (z1 - z2) * (y - b_y1) / (b_y1 - y2) + z1;

            /* 'find_intercepts_3d_4dof_dfki:57' i_num = i_num + 1; */
            *i_num = 2.0;
        } else if (floor(b_y1 / depth_cell_size) == ceil(y2 / depth_cell_size) + 1.0)
        {
            /* 'find_intercepts_3d_4dof_dfki:58' elseif floor(y1/depth_cell_size) == ceil(y2/depth_cell_size)+1 */
            /*  2 depth_cell intercepts */
            /* 'find_intercepts_3d_4dof_dfki:60' y = floor(y1/depth_cell_size)*depth_cell_size; */
            y = floor(b_y1 / depth_cell_size) * depth_cell_size;

            /* 'find_intercepts_3d_4dof_dfki:61' x = (x1-x2)*(y-y1)/(y1-y2) + x1; */
            /* 'find_intercepts_3d_4dof_dfki:62' z = (z1-z2)*(y-y1)/(y1-y2) + z1; */
            /* 'find_intercepts_3d_4dof_dfki:63' i(i_num,1) = x; */
            i[0] = (x1 - x2) * (y - b_y1) / (b_y1 - y2) + x1;

            /* 'find_intercepts_3d_4dof_dfki:64' i(i_num,2) = y; */
            i[7] = y;

            /* 'find_intercepts_3d_4dof_dfki:65' i(i_num,3) = z; */
            i[14] = (z1 - z2) * (y - b_y1) / (b_y1 - y2) + z1;

            /* 'find_intercepts_3d_4dof_dfki:66' i_num = i_num + 1; */
            /* 'find_intercepts_3d_4dof_dfki:68' y = ceil(y2/depth_cell_size)*depth_cell_size; */
            y = ceil(y2 / depth_cell_size) * depth_cell_size;

            /* 'find_intercepts_3d_4dof_dfki:69' x = (x1-x2)*(y-y1)/(y1-y2) + x1; */
            /* 'find_intercepts_3d_4dof_dfki:70' z = (z1-z2)*(y-y1)/(y1-y2) + z1; */
            /* 'find_intercepts_3d_4dof_dfki:71' i(i_num,1) = x; */
            i[1] = (x1 - x2) * (y - b_y1) / (b_y1 - y2) + x1;

            /* 'find_intercepts_3d_4dof_dfki:72' i(i_num,2) = y; */
            i[8] = y;

            /* 'find_intercepts_3d_4dof_dfki:73' i(i_num,3) = z; */
            i[15] = (z1 - z2) * (y - b_y1) / (b_y1 - y2) + z1;

            /* 'find_intercepts_3d_4dof_dfki:74' i_num = i_num + 1; */
            *i_num = 3.0;
        } else {
            /* 'find_intercepts_3d_4dof_dfki:76' else */
            /*  no depth_cell_intercepts */
        }
    }

    /*  hori_cell intercept (can only have one if depth_cell_size is same as the measurement cell) */
    /* 'find_intercepts_3d_4dof_dfki:83' if floor(x2/hori_res) == ceil(x1/hori_res) */
    if (floor(x2 / hori_res) == ceil(x1 / hori_res)) {
        /* 'find_intercepts_3d_4dof_dfki:84' x = floor(x2/hori_res)*hori_res; */
        y = floor(x2 / hori_res) * hori_res;

        /* 'find_intercepts_3d_4dof_dfki:85' y = (y1-y2)*(x-x1)/(x1-x2) + y1; */
        /* 'find_intercepts_3d_4dof_dfki:86' z = (z1-z2)*(x-x1)/(x1-x2) + z1; */
        /* 'find_intercepts_3d_4dof_dfki:87' i(i_num,1) = x; */
        i[(int32_T)*i_num - 1] = y;

        /* 'find_intercepts_3d_4dof_dfki:88' i(i_num,2) = y; */
        i[(int32_T)*i_num + 6] = (b_y1 - y2) * (y - x1) / (x1 - x2) + b_y1;

        /* 'find_intercepts_3d_4dof_dfki:89' i(i_num,3) = z; */
        i[(int32_T)*i_num + 13] = (z1 - z2) * (y - x1) / (x1 - x2) + z1;

        /* 'find_intercepts_3d_4dof_dfki:90' i_num = i_num + 1; */
        (*i_num)++;

        /*  1 hori_cell intercept in x */
    }

    /* 'find_intercepts_3d_4dof_dfki:95' if floor(z2/hori_res) == ceil(z1/hori_res) */
    if (floor(z2 / hori_res) == ceil(z1 / hori_res)) {
        /* 'find_intercepts_3d_4dof_dfki:96' z = floor(z2/hori_res)*hori_res; */
        y = floor(z2 / hori_res) * hori_res;

        /* 'find_intercepts_3d_4dof_dfki:97' x = (x1-x2)*(z-z1)/(z1-z2) + x1; */
        /* 'find_intercepts_3d_4dof_dfki:98' y = (y1-y2)*(z-z1)/(z1-z2) + y1; */
        /* 'find_intercepts_3d_4dof_dfki:99' i(i_num,1) = x; */
        i[(int32_T)*i_num - 1] = (x1 - x2) * (y - z1) / (z1 - z2) + x1;

        /* 'find_intercepts_3d_4dof_dfki:100' i(i_num,2) = y; */
        i[(int32_T)*i_num + 6] = (b_y1 - y2) * (y - z1) / (z1 - z2) + b_y1;

        /* 'find_intercepts_3d_4dof_dfki:101' i(i_num,3) = z; */
        i[(int32_T)*i_num + 13] = y;

        /* 'find_intercepts_3d_4dof_dfki:102' i_num = i_num + 1; */
        (*i_num)++;

        /*  1 hori_cell intercept in z */
    }

    /* 'find_intercepts_3d_4dof_dfki:106' if  ceil(z2/hori_res) == floor(z1/hori_res) */
    if (ceil(z2 / hori_res) == floor(z1 / hori_res)) {
        /* 'find_intercepts_3d_4dof_dfki:107' z = ceil(z2/hori_res)*hori_res; */
        y = ceil(z2 / hori_res) * hori_res;

        /* 'find_intercepts_3d_4dof_dfki:108' x = (x1-x2)*(z-z1)/(z1-z2) + x1; */
        /* 'find_intercepts_3d_4dof_dfki:109' y = (y1-y2)*(z-z1)/(z1-z2) + y1; */
        /* 'find_intercepts_3d_4dof_dfki:110' i(i_num,1) = x; */
        i[(int32_T)*i_num - 1] = (x1 - x2) * (y - z1) / (z1 - z2) + x1;

        /* 'find_intercepts_3d_4dof_dfki:111' i(i_num,2) = y; */
        i[(int32_T)*i_num + 6] = (b_y1 - y2) * (y - z1) / (z1 - z2) + b_y1;

        /* 'find_intercepts_3d_4dof_dfki:112' i(i_num,3) = z; */
        i[(int32_T)*i_num + 13] = y;

        /* 'find_intercepts_3d_4dof_dfki:113' i_num = i_num + 1; */
        (*i_num)++;

        /*  1 hori_cell intercept in z */
    }

    /*  end */
    /*  if i_num == 1 */
    /*     i = [ ]; */
    /*  end */
}

/* End of code generation (find_intercepts_3d_4dof_dfki.cpp) */

/*
               * function w = find_weighting_trilin(x0,y0,z0,xe,ye,ze,l1,l2)
               */
void ADCP_measurement_model::find_weighting_trilin(real_T x0, real_T b_y0, real_T z0, real_T xe, real_T
                                                   ye, real_T ze, real_T l1, real_T l2, real_T w[8])
{
    real_T Wv8l2l1;
    real_T y;
    real_T b_y;
    real_T c_y;
    real_T d_y;
    real_T e_y;
    real_T Wv1l2l1;
    real_T Wv2l2l1;
    real_T Wv3l2l1;
    real_T Wv4l2l1;
    real_T Wv5l2l1;
    real_T Wv6l2l1;
    real_T Wv7l2l1;

    /*  x1 is the first cut of the triangle */
    /*  x2 is the second cut of the triangle */
    /*  overhang is the amount which overhangs to the next measurement cell */
    /*  depth_cell_size is the size of the measurement cell */
    /*  w is the area of the triangle (normalized from 0 to 1) for the weighting */
    /* 'find_weighting_trilin:8' if l2 < l1 */
    /* 'find_weighting_trilin:12' if l1 <= 1/2 && l2 <= 1/2 */
    if ((l1 <= 0.5) && (l2 <= 0.5)) {
        /* 'find_weighting_trilin:14' Wv1l2l1 = 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) + (x0 - xe)*(y0 - ye)*(ze - 1))*l1^4 + (4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) + 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l1^3 + 2*(xe - 1)*(ye - 1)*(ze - 1)*l1^2 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- ((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - (x0 - xe)*(y0 - ye)*(ze - 1))*l2^4 + (- 4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) - 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l2^3 - 2*(xe - 1)*(ye - 1)*(ze - 1)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv1l2l1 = ((((((0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 + (((x0
                                                                                 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) + (x0 - xe) *
                                                                               (b_y0 - ye) * (ze - 1.0)) * y) + (1.3333333333333333 * ((x0 - xe) * (ye -
                                                                                                                                                    1.0) + (b_y0 - ye) * (xe - 1.0)) * (ze - 1.0) + 1.3333333333333333 * (z0 -
                                                                                                                                                                                                                          ze) * (xe - 1.0) * (ye - 1.0)) * b_y) + 2.0 * (xe - 1.0) * (ye - 1.0) *
                      (ze - 1.0) * (l1 * l1)) - 0.8 * (x0 - xe) * (b_y0 - ye) * (z0
                                                                                 - ze) * c_y) + (-((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0
                                                                                                                                                         - ze) - (x0 - xe) * (b_y0 - ye) * (ze - 1.0)) * d_y) +
                   (-1.3333333333333333 * ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) *
                                           (xe - 1.0)) * (ze - 1.0) - 1.3333333333333333 * (z0 - ze) * (xe - 1.0) *
                    (ye - 1.0)) * e_y) - 2.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0) *
                (l2 * l2);

        /* 'find_weighting_trilin:15' Wv2l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1))*(x0 - xe) - xe*(y0 - ye)*(z0 - ze))*l1^4 + (- 1/3*xe*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) - 4/3*(x0 - xe)*(ye - 1)*(ze - 1))*l1^3 - 2*xe*(ye - 1)*(ze - 1)*l1^2 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1))*(x0 - xe) + xe*(y0 - ye)*(z0 - ze))*l2^4 + (1/3*xe*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) + 4/3*(x0 - xe)*(ye - 1)*(ze - 1))*l2^3 + 2*xe*(ye - 1)*(ze - 1)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv2l2l1 = ((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                        (-0.25 * (4.0 * (b_y0 - ye) * (ze - 1.0) + 4.0 * (z0 - ze) *
                                  (ye - 1.0)) * (x0 - xe) - xe * (b_y0 - ye) * (z0 - ze)) *
                        y) + (-0.33333333333333331 * xe * (4.0 * (b_y0 - ye) * (ze -
                                                                                1.0) + 4.0 * (z0 - ze) * (ye - 1.0)) - 1.3333333333333333 * (x0 - xe) *
                              (ye - 1.0) * (ze - 1.0)) * b_y) - 2.0 * xe * (ye - 1.0)
                      * (ze - 1.0) * (l1 * l1)) + 0.8 * (x0 - xe) * (b_y0 - ye) *
                     (z0 - ze) * c_y) + (0.25 * (4.0 * (b_y0 - ye) * (ze - 1.0) +
                                                 4.0 * (z0 - ze) * (ye - 1.0)) * (x0 - xe) + xe * (b_y0 - ye) * (z0 - ze)) *
                    d_y) + (0.33333333333333331 * xe * (4.0 * (b_y0 - ye) * (ze -
                                                                             1.0) + 4.0 * (z0 - ze) * (ye - 1.0)) + 1.3333333333333333 * (x0 - xe) *
                            (ye - 1.0) * (ze - 1.0)) * e_y) + 2.0 * xe * (ye - 1.0) *
                (ze - 1.0) * (l2 * l2);

        /* 'find_weighting_trilin:16' Wv3l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1))*(y0 - ye) - ye*(x0 - xe)*(z0 - ze))*l1^4 + (- 1/3*ye*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) - 4/3*(y0 - ye)*(xe - 1)*(ze - 1))*l1^3 - 2*ye*(xe - 1)*(ze - 1)*l1^2 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1))*(y0 - ye) + ye*(x0 - xe)*(z0 - ze))*l2^4 + (1/3*ye*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) + 4/3*(y0 - ye)*(xe - 1)*(ze - 1))*l2^3 + 2*ye*(xe - 1)*(ze - 1)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv3l2l1 = ((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                        (-0.25 * (4.0 * (x0 - xe) * (ze - 1.0) + 4.0 * (z0 - ze) *
                                  (xe - 1.0)) * (b_y0 - ye) - ye * (x0 - xe) * (z0 - ze)) * y) +
                       (-0.33333333333333331 * ye * (4.0 * (x0 - xe) * (ze - 1.0) +
                                                     4.0 * (z0 - ze) * (xe - 1.0)) - 1.3333333333333333 * (b_y0 - ye) * (xe -
                                                                                                                         1.0) * (ze - 1.0)) * b_y) - 2.0 * ye * (xe - 1.0) * (ze - 1.0) * (l1 * l1))
                     + 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y) + (0.25 *
                                                                           (4.0 * (x0 - xe) * (ze - 1.0) + 4.0 * (z0 - ze) * (xe - 1.0)) * (b_y0 - ye)
                                                                           + ye * (x0 - xe) * (z0 - ze)) * d_y) + (0.33333333333333331 * ye * (4.0 *
                                                                                                                                               (x0 - xe) * (ze - 1.0) + 4.0 * (z0 - ze) * (xe - 1.0)) +
                                                                                                                   1.3333333333333333 * (b_y0 - ye) * (xe - 1.0) * (ze - 1.0)) *
                   e_y) + 2.0 * ye * (xe - 1.0) * (ze - 1.0) * (l2 * l2);

        /* 'find_weighting_trilin:17' Wv4l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1))*(z0 - ze) - ze*(x0 - xe)*(y0 - ye))*l1^4 + (- 1/3*ze*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) - 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l1^3 - 2*ze*(xe - 1)*(ye - 1)*l1^2 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1))*(z0 - ze) + ze*(x0 - xe)*(y0 - ye))*l2^4 + (1/3*ze*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) + 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l2^3 + 2*ze*(xe - 1)*(ye - 1)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv4l2l1 = ((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                        (-0.25 * (4.0 * (x0 - xe) * (ye - 1.0) + 4.0 * (b_y0 - ye) *
                                  (xe - 1.0)) * (z0 - ze) - ze * (x0 - xe) * (b_y0 - ye)) *
                        y) + (-0.33333333333333331 * ze * (4.0 * (x0 - xe) * (ye -
                                                                              1.0) + 4.0 * (b_y0 - ye) * (xe - 1.0)) - 1.3333333333333333 * (z0 - ze) *
                              (xe - 1.0) * (ye - 1.0)) * b_y) - 2.0 * ze * (xe - 1.0)
                      * (ye - 1.0) * (l1 * l1)) + 0.8 * (x0 - xe) * (b_y0 - ye) *
                     (z0 - ze) * c_y) + (0.25 * (4.0 * (x0 - xe) * (ye - 1.0) + 4.0 *
                                                 (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) + ze * (x0 - xe) * (b_y0
                                                                                                           - ye)) * d_y) + (0.33333333333333331 * ze * (4.0 * (x0 - xe) * (ye - 1.0)
                                                                                                                                                        + 4.0 * (b_y0 - ye) * (xe - 1.0)) + 1.3333333333333333 * (z0 - ze) * (xe -
                                                                                                                                                                                                                              1.0) * (ye - 1.0)) * e_y) + 2.0 * ze * (xe - 1.0) * (ye - 1.0) * (l2 * l2);

        /* 'find_weighting_trilin:18' Wv5l2l1 = 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l1^5 + (1/4*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye))*(z0 - ze) + 1/4*ze*(4*y0 - 4*ye)*(x0 - xe))*l1^4 + (1/3*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye)) + 1/3*xe*(4*ye - 4)*(z0 - ze))*l1^3 + 1/2*xe*ze*(4*ye - 4)*l1^2 - 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l2^5 + (- 1/4*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye))*(z0 - ze) - 1/4*ze*(4*y0 - 4*ye)*(x0 - xe))*l2^4 + (- 1/3*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye)) - 1/3*xe*(4*ye - 4)*(z0 - ze))*l2^3 - 1/2*xe*ze*(4*ye - 4)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv5l2l1 = ((((((0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 - ze) *
                        Wv8l2l1 + (0.25 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * (4.0 *
                                                                                b_y0 - 4.0 * ye)) * (z0 - ze) + 0.25 * ze * (4.0 * b_y0 - 4.0 * ye) * (x0
                                                                                                                                                       - xe)) * y) + (0.33333333333333331 * ze * ((4.0 * ye - 4.0) * (x0 - xe) +
                                                                                                                                                                                                  xe * (4.0 * b_y0 - 4.0 * ye)) + 0.33333333333333331 * xe * (4.0 * ye - 4.0)
                                                                                                                                                                      * (z0 - ze)) * b_y) + 0.5 * xe * ze * (4.0 * ye - 4.0) *
                      (l1 * l1)) - 0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 -
                                                                                ze) * c_y) + (-0.25 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * (4.0 * b_y0 -
                                                                                                                                            4.0 * ye)) * (z0 - ze) - 0.25 * ze * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe)) *
                    d_y) + (-0.33333333333333331 * ze * ((4.0 * ye - 4.0) * (x0 - xe)
                                                         + xe * (4.0 * b_y0 - 4.0 * ye)) - 0.33333333333333331 * xe * (4.0 * ye -
                                                                                                                       4.0) * (z0 - ze)) * e_y) - 0.5 * xe * ze * (4.0 * ye - 4.0) * (l2 * l2);

        /* 'find_weighting_trilin:19' Wv6l2l1 = 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe))*(z0 - ze) + 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l1^4 + (1/3*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe)) + 1/3*ye*(4*xe - 4)*(z0 - ze))*l1^3 + 1/2*ye*ze*(4*xe - 4)*l1^2 - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe))*(z0 - ze) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l2^4 + (- 1/3*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe)) - 1/3*ye*(4*xe - 4)*(z0 - ze))*l2^3 - 1/2*ye*ze*(4*xe - 4)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv6l2l1 = ((((((0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                        Wv8l2l1 + (0.25 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye *
                                           (4.0 * x0 - 4.0 * xe)) * (z0 - ze) + 0.25 * ze * (4.0 * x0 - 4.0 * xe) *
                                   (b_y0 - ye)) * y) + (0.33333333333333331 * ze * ((4.0 * xe - 4.0) * (b_y0
                                                                                                        - ye) + ye * (4.0 * x0 - 4.0 * xe)) + 0.33333333333333331 * ye * (4.0 * xe
                                                                                                                                                                          - 4.0) * (z0 - ze)) * b_y) + 0.5 * ye * ze * (4.0 * xe - 4.0) * (l1 * l1))
                     - 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) * c_y)
                    + (-0.25 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * (4.0 * x0 -
                                                                       4.0 * xe)) * (z0 - ze) - 0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye)) *
                    d_y) + (-0.33333333333333331 * ze * ((4.0 * xe - 4.0) * (b_y0 -
                                                                             ye) + ye * (4.0 * x0 - 4.0 * xe)) - 0.33333333333333331 * ye * (4.0 * xe -
                                                                                                                                             4.0) * (z0 - ze)) * e_y) - 0.5 * ye * ze * (4.0 * xe - 4.0) * (l2 * l2);

        /* 'find_weighting_trilin:20' Wv7l2l1 = 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l1^5 + (1/4*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze))*(y0 - ye) + 1/4*ye*(4*z0 - 4*ze)*(x0 - xe))*l1^4 + (1/3*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze)) + 1/3*xe*(4*ze - 4)*(y0 - ye))*l1^3 + 1/2*xe*ye*(4*ze - 4)*l1^2 - 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l2^5 + (- 1/4*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze))*(y0 - ye) - 1/4*ye*(4*z0 - 4*ze)*(x0 - xe))*l2^4 + (- 1/3*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze)) - 1/3*xe*(4*ze - 4)*(y0 - ye))*l2^3 - 1/2*xe*ye*(4*ze - 4)*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv7l2l1 = ((((((0.2 * (4.0 * z0 - 4.0 * ze) * (x0 - xe) * (b_y0 - ye) *
                        Wv8l2l1 + (0.25 * ((4.0 * ze - 4.0) * (x0 - xe) + xe * (4.0 *
                                                                                z0 - 4.0 * ze)) * (b_y0 - ye) + 0.25 * ye * (4.0 * z0 - 4.0 * ze) * (x0 -
                                                                                                                                                     xe)) * y) + (0.33333333333333331 * ye * ((4.0 * ze - 4.0) * (x0 - xe) + xe
                                                                                                                                                                                              * (4.0 * z0 - 4.0 * ze)) + 0.33333333333333331 * xe * (4.0 * ze - 4.0) *
                                                                                                                                                                  (b_y0 - ye)) * b_y) + 0.5 * xe * ye * (4.0 * ze - 4.0) * (l1 *
                                                                                                                                                                                                                            l1)) - 0.2 * (4.0 * z0 - 4.0 * ze) * (x0 - xe) * (b_y0 - ye) * c_y) +
                    (-0.25 * ((4.0 * ze - 4.0) * (x0 - xe) + xe * (4.0 * z0 - 4.0 *
                                                                   ze)) * (b_y0 - ye) - 0.25 * ye * (4.0 * z0 - 4.0 * ze) * (x0 - xe)) * d_y)
                   + (-0.33333333333333331 * ye * ((4.0 * ze - 4.0) * (x0 - xe) + xe
                                                   * (4.0 * z0 - 4.0 * ze)) - 0.33333333333333331 * xe * (4.0 * ze - 4.0) *
                      (b_y0 - ye)) * e_y) - 0.5 * xe * ye * (4.0 * ze - 4.0) * (l2 *
                                                                                l2);

        /* 'find_weighting_trilin:21' Wv8l2l1 = - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l1^4 + (- 1/3*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) - 4/3*xe*ye*(z0 - ze))*l1^3 - 2*xe*ye*ze*l1^2 + 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) + 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l2^4 + (1/3*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) + 4/3*xe*ye*(z0 - ze))*l2^3 + 2*xe*ye*ze*l2^2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv8l2l1 = ((((((-0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                        Wv8l2l1 + (-0.25 * (z0 - ze) * (4.0 * xe * (b_y0 - ye) + ye *
                                                        (4.0 * x0 - 4.0 * xe)) - 0.25 * ze * (4.0 * x0 - 4.0 * xe)
                                   * (b_y0 - ye)) * y) + (-0.33333333333333331 * ze * (4.0 * xe * (b_y0 - ye)
                                                                                       + ye * (4.0 * x0 - 4.0 * xe)) - 1.3333333333333333 * xe * ye * (z0 - ze)) *
                       b_y) - 2.0 * xe * ye * ze * (l1 * l1)) + 0.2 * (4.0 * x0 -
                                                                       4.0 * xe) * (b_y0 - ye) * (z0 - ze) * c_y) + (0.25 * (z0 - ze) * (4.0 * xe
                                                                                                                                         * (b_y0 - ye) + ye * (4.0 * x0 - 4.0 * xe)) + 0.25 * ze * (4.0 * x0 - 4.0 *
                                                                                                                                                                                                    xe) * (b_y0 - ye)) * d_y) + (0.33333333333333331 * ze * (4.0 * xe * (b_y0
                                                                                                                                                                                                                                                                         - ye) + ye * (4.0 * x0 - 4.0 * xe)) + 1.3333333333333333 * xe * ye * (z0 -
                                                                                                                                                                                                                                                                                                                                               ze)) * e_y) + 2.0 * xe * ye * ze * (l2 * l2);
    } else if ((l1 >= 0.5) && (l2 >= 0.5)) {
        /* 'find_weighting_trilin:23' elseif l1 >= 1/2 && l2 >= 1/2 */
        /* 'find_weighting_trilin:25' Wv1l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + ((x0 - xe)*(y0 - ye)*(z0 - ze) - ((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - (x0 - xe)*(y0 - ye)*(ze - 1))*l1^4 + (4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - 4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) + 4/3*(x0 - xe)*(y0 - ye)*(ze - 1) - 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l1^3 + (2*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) + 2*(z0 - ze)*(xe - 1)*(ye - 1) - 2*(xe - 1)*(ye - 1)*(ze - 1))*l1^2 + 4*(xe - 1)*(ye - 1)*(ze - 1)*l1 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - (x0 - xe)*(y0 - ye)*(z0 - ze) + (x0 - xe)*(y0 - ye)*(ze - 1))*l2^4 + (4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) - 4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - 4/3*(x0 - xe)*(y0 - ye)*(ze - 1) + 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l2^3 + (2*(xe - 1)*(ye - 1)*(ze - 1) - 2*(z0 - ze)*(xe - 1)*(ye - 1) - 2*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1))*l2^2 - 4*(xe - 1)*(ye - 1)*(ze - 1)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv1l2l1 = ((((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                          (((x0 - xe) * (b_y0 - ye) * (z0 - ze) - ((x0 - xe) * (ye -
                                                                                1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0 - ze)) - (x0 - xe) * (b_y0 - ye) *
                           (ze - 1.0)) * y) + (((1.3333333333333333 * ((x0 - xe) *
                                                                       (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) - 1.3333333333333333 *
                                                 ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) *
                                                 (ze - 1.0)) + 1.3333333333333333 * (x0 - xe) * (b_y0 - ye) * (ze - 1.0)) -
                                               1.3333333333333333 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) * b_y) + ((2.0 *
                                                                                                                    ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (ze
                                                                                                                                                                           - 1.0) + 2.0 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) - 2.0 * (xe - 1.0) *
                                                                                                                   (ye - 1.0) * (ze - 1.0)) * (l1 * l1)) + 4.0 * (xe - 1.0) * (ye - 1.0) *
                       (ze - 1.0) * l1) + 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) *
                      c_y) + ((((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) *
                               (z0 - ze) - (x0 - xe) * (b_y0 - ye) * (z0 - ze)) + (x0 - xe) *
                              (b_y0 - ye) * (ze - 1.0)) * d_y) +
                    (((1.3333333333333333 * ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) *
                                             (xe - 1.0)) * (ze - 1.0) - 1.3333333333333333 * ((x0 - xe) *
                                                                                              (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0 - ze)) -
                      1.3333333333333333 * (x0 - xe) * (b_y0 - ye) * (ze - 1.0)) +
                     1.3333333333333333 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) * e_y)
                   + ((2.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0) - 2.0 * (z0 - ze) *
                       (xe - 1.0) * (ye - 1.0)) - 2.0 * ((x0 - xe) * (ye - 1.0) +
                                                         (b_y0 - ye) * (xe - 1.0)) * (ze - 1.0)) * (l2 * l2)) - 4.0 * (xe - 1.0) *
                (ye - 1.0) * (ze - 1.0) * l2;

        /* 'find_weighting_trilin:26' Wv2l2l1 = 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)) + xe*(y0 - ye)*(z0 - ze))*l1^4 + (1/3*xe*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)) - 1/3*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)))*l1^3 + (- 1/2*xe*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) - 2*(x0 - xe)*(ye - 1)*(ze - 1))*l1^2 - 4*xe*(ye - 1)*(ze - 1)*l1 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)) - xe*(y0 - ye)*(z0 - ze))*l2^4 + (1/3*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) - 1/3*xe*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)))*l2^3 + (1/2*xe*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) + 2*(x0 - xe)*(ye - 1)*(ze - 1))*l2^2 + 4*xe*(ye - 1)*(ze - 1)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv2l2l1 = ((((((((0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                          (0.25 * (x0 - xe) * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0
                                                * (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0 - ze) * (ye - 1.0)) + xe * (b_y0 -
                                                                                                                   ye) * (z0 - ze)) * y) + (0.33333333333333331 * xe * ((4.0 * (b_y0 - ye) *
                                                                                                                                                                         (ze - 1.0) - 4.0 * (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0
                                                                                                                                                                                                                              - ze) * (ye - 1.0)) - 0.33333333333333331 * (x0 - xe) * ((4.0 * (b_y0 - ye)
                                                                                                                                                                                                                                                                                        * (ze - 1.0) - 4.0 * (ye - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (ye -
                                                                                                                                                                                                                                                                                                                                                           1.0))) * b_y) + (-0.5 * xe * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0 * (ye
                                                                                                                                                                                                                                                                                                                                                                                                                                  - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (ye - 1.0)) - 2.0 * (x0 - xe) *
                                                                                                                                                                                                                                                                                                                                                                            (ye - 1.0) * (ze - 1.0)) * (l1 * l1)) - 4.0 * xe * (ye -
                                                                                                                                                                                                                                                                                                                                                                                                                                1.0) * (ze - 1.0) * l1) - 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y)
                     + (-0.25 * (x0 - xe) * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0 *
                                              (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0 - ze) * (ye - 1.0)) -
                        xe * (b_y0 - ye) * (z0 - ze)) * d_y) + (0.33333333333333331 *
                                                                (x0 - xe) * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0 * (ye - 1.0)
                                                                              * (ze - 1.0)) + 4.0 * (z0 - ze) * (ye - 1.0)) - 0.33333333333333331 * xe *
                                                                ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0 * (b_y0 - ye) * (z0 - ze))
                                                                 + 4.0 * (z0 - ze) * (ye - 1.0))) * e_y) + (0.5 * xe * ((4.0 *
                                                                                                                         (b_y0 - ye) * (ze - 1.0) - 4.0 * (ye - 1.0) * (ze - 1.0)) +
                                                                                                                        4.0 * (z0 - ze) * (ye - 1.0)) + 2.0 * (x0 - xe) * (ye - 1.0) * (ze - 1.0))
                   * (l2 * l2)) + 4.0 * xe * (ye - 1.0) * (ze - 1.0) * l2;

        /* 'find_weighting_trilin:27' Wv3l2l1 = 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)) + ye*(x0 - xe)*(z0 - ze))*l1^4 + (1/3*ye*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)) - 1/3*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)))*l1^3 + (- 1/2*ye*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) - 2*(y0 - ye)*(xe - 1)*(ze - 1))*l1^2 - 4*ye*(xe - 1)*(ze - 1)*l1 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)) - ye*(x0 - xe)*(z0 - ze))*l2^4 + (1/3*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) - 1/3*ye*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)))*l2^3 + (1/2*ye*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) + 2*(y0 - ye)*(xe - 1)*(ze - 1))*l2^2 + 4*ye*(xe - 1)*(ze - 1)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv3l2l1 = ((((((((0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                          (0.25 * (b_y0 - ye) * ((4.0 * (x0 - xe) * (ze - 1.0) - 4.0
                                                  * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) * (xe - 1.0)) + ye * (x0 - xe) *
                           (z0 - ze)) * y) + (0.33333333333333331 * ye * ((4.0 * (x0
                                                                                  - xe) * (ze - 1.0) - 4.0 * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) * (xe
                                                                                                                                                         - 1.0)) - 0.33333333333333331 * (b_y0 - ye) * ((4.0 * (x0 - xe) * (ze -
                                                                                                                                                                                                                            1.0) - 4.0 * (xe - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (xe - 1.0))) *
                         b_y) + (-0.5 * ye * ((4.0 * (x0 - xe) * (ze - 1.0) - 4.0 *
                                               (xe - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (xe - 1.0))
                                 - 2.0 * (b_y0 - ye) * (xe - 1.0) * (ze - 1.0)) * (l1 * l1)) - 4.0 * ye *
                       (xe - 1.0) * (ze - 1.0) * l1) - 0.8 * (x0 - xe) * (b_y0 - ye)
                      * (z0 - ze) * c_y) + (-0.25 * (b_y0 - ye) * ((4.0 * (x0 - xe) *
                                                                    (ze - 1.0) - 4.0 * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) *
                                                                   (xe - 1.0)) - ye * (x0 - xe) * (z0 - ze)) * d_y) +
                    (0.33333333333333331 * (b_y0 - ye) * ((4.0 * (x0 - xe) * (ze -
                                                                              1.0) - 4.0 * (xe - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (xe - 1.0)) -
                     0.33333333333333331 * ye * ((4.0 * (x0 - xe) * (ze - 1.0) - 4.0
                                                  * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) * (xe - 1.0))) * e_y) + (0.5 *
                                                                                                                      ye * ((4.0 * (x0 - xe) * (ze - 1.0) - 4.0 * (xe - 1.0) * (ze -
                                                                                                                                                                                1.0)) + 4.0 * (z0 - ze) * (xe - 1.0)) + 2.0 * (b_y0 - ye) * (xe - 1.0) *
                                                                                                                      (ze - 1.0)) * (l2 * l2)) + 4.0 * ye * (xe - 1.0) * (ze - 1.0) *
                l2;

        /* 'find_weighting_trilin:28' Wv4l2l1 = 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)) + ze*(x0 - xe)*(y0 - ye))*l1^4 + (1/3*ze*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)) - 1/3*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)))*l1^3 + (- 1/2*ze*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) - 2*(z0 - ze)*(xe - 1)*(ye - 1))*l1^2 - 4*ze*(xe - 1)*(ye - 1)*l1 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)) - ze*(x0 - xe)*(y0 - ye))*l2^4 + (1/3*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) - 1/3*ze*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)))*l2^3 + (1/2*ze*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) + 2*(z0 - ze)*(xe - 1)*(ye - 1))*l2^2 + 4*ze*(xe - 1)*(ye - 1)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv4l2l1 = ((((((((0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1 +
                          (0.25 * (z0 - ze) * ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 *
                                                (x0 - xe) * (b_y0 - ye)) + 4.0 * (b_y0 - ye) * (xe -
                                                                                                1.0)) + ze * (x0 - xe) * (b_y0 - ye)) * y) + (0.33333333333333331 * ze *
                                                                                                                                              ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (x0 - xe) * (b_y0 -
                                                                                                                                                                                                  ye)) + 4.0 * (b_y0 - ye) * (xe - 1.0)) - 0.33333333333333331 * (z0 - ze) *
                                                                                                                                              ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (xe - 1.0) * (ye -
                                                                                                                                                                                                   1.0)) + 4.0 * (b_y0 - ye) * (xe - 1.0))) * b_y) + (-0.5 * ze * ((4.0 * (x0
                                                                                                                                                                                                                                                                           - xe) * (ye - 1.0) - 4.0 * (xe - 1.0) * (ye - 1.0)) + 4.0 * (b_y0 - ye) *
                                                                                                                                                                                                                                                                   (xe - 1.0)) - 2.0 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) *
                        (l1 * l1)) - 4.0 * ze * (xe - 1.0) * (ye - 1.0) * l1) - 0.8 *
                      (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y) + (-0.25 * (z0 - ze)
                                                                    * ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (x0 - xe) * (b_y0 - ye)) + 4.0 *
                                                                       (b_y0 - ye) * (xe - 1.0)) - ze * (x0 - xe) * (b_y0 - ye)) * d_y) +
                    (0.33333333333333331 * (z0 - ze) * ((4.0 * (x0 - xe) * (ye - 1.0)
                                                         - 4.0 * (xe - 1.0) * (ye - 1.0)) + 4.0 * (b_y0 - ye) * (xe - 1.0)) -
                     0.33333333333333331 * ze * ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0
                                                  * (x0 - xe) * (b_y0 - ye)) + 4.0 * (b_y0 - ye) * (xe - 1.0))) * e_y) +
                   (0.5 * ze * ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (xe - 1.0) *
                                 (ye - 1.0)) + 4.0 * (b_y0 - ye) * (xe - 1.0)) + 2.0 * (z0 - ze) * (xe -
                                                                                                    1.0) * (ye - 1.0)) * (l2 * l2)) + 4.0 * ze * (xe - 1.0) * (ye - 1.0) * l2;

        /* 'find_weighting_trilin:29' Wv5l2l1 = - 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l1^5 + (1/4*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye))*(z0 - ze) - 1/4*ze*(4*y0 - 4*ye)*(x0 - xe))*l1^4 + (1/3*ze*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye)) + 1/3*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4))*(z0 - ze))*l1^3 + (1/2*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4)) + 1/2*xe*(4*ye - 4)*(z0 - ze))*l1^2 + xe*ze*(4*ye - 4)*l1 + 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l2^5 + (1/4*ze*(4*y0 - 4*ye)*(x0 - xe) - 1/4*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye))*(z0 - ze))*l2^4 + (- 1/3*ze*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye)) - 1/3*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4))*(z0 - ze))*l2^3 + (- 1/2*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4)) - 1/2*xe*(4*ye - 4)*(z0 - ze))*l2^2 - xe*ze*(4*ye - 4)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv5l2l1 = ((((((((-0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 - ze) *
                          Wv8l2l1 + (0.25 * ((x0 - xe) * ((4.0 * b_y0 - 8.0 * ye) +
                                                          4.0) - xe * (4.0 * b_y0 - 4.0 * ye)) * (z0 - ze) - 0.25 * ze * (4.0 * b_y0
                                                                                                                          - 4.0 * ye) * (x0 - xe)) * y) + (0.33333333333333331 * ze * ((x0 - xe) *
                                                                                                                                                                                       ((4.0 * b_y0 - 8.0 * ye) + 4.0) - xe * (4.0 * b_y0 - 4.0 *
                                                                                                                                                                                                                               ye)) + 0.33333333333333331 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * ((4.0 *
                                                                                                                                                                                                                                                                                                   b_y0 - 8.0 * ye) + 4.0)) * (z0 - ze)) * b_y) + (0.5 * ze * ((4.0 * ye -
                                                                                                                                                                                                                                                                                                                                                                4.0) * (x0 - xe) + xe * ((4.0 * b_y0 - 8.0 * ye) + 4.0)) + 0.5 * xe * (4.0
                                                                                                                                                                                                                                                                                                                                                                                                                                       * ye - 4.0) * (z0 - ze)) * (l1 * l1)) + xe * ze * (4.0 * ye - 4.0) * l1) +
                      0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 - ze) * c_y) +
                     (0.25 * ze * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) - 0.25 * ((x0
                                                                                 - xe) * ((4.0 * b_y0 - 8.0 * ye) + 4.0) - xe * (4.0 * b_y0 - 4.0 * ye)) *
                      (z0 - ze)) * d_y) + (-0.33333333333333331 * ze * ((x0 - xe) *
                                                                        ((4.0 * b_y0 - 8.0 * ye) + 4.0) - xe * (4.0 * b_y0 - 4.0 * ye))
                                           - 0.33333333333333331 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * ((4.0 * b_y0
                                                                                                          - 8.0 * ye) + 4.0)) * (z0 - ze)) * e_y) + (-0.5 * ze * ((4.0 * ye - 4.0) *
                                                                                                                                                                  (x0 - xe) + xe * ((4.0 * b_y0 - 8.0 * ye) + 4.0)) - 0.5 * xe *
                                                                                                                                                     (4.0 * ye - 4.0) * (z0 - ze)) * (l2 * l2)) - xe * ze * (4.0 * ye
                                                                                                                                                                                                             - 4.0) * l2;

        /* 'find_weighting_trilin:30' Wv6l2l1 = - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe))*(z0 - ze) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l1^4 + (1/3*ze*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe)) + 1/3*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4))*(z0 - ze))*l1^3 + (1/2*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4)) + 1/2*ye*(4*xe - 4)*(z0 - ze))*l1^2 + ye*ze*(4*xe - 4)*l1 + 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*ze*(4*x0 - 4*xe)*(y0 - ye) - 1/4*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe))*(z0 - ze))*l2^4 + (- 1/3*ze*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe)) - 1/3*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4))*(z0 - ze))*l2^3 + (- 1/2*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4)) - 1/2*ye*(4*xe - 4)*(z0 - ze))*l2^2 - ye*ze*(4*xe - 4)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv6l2l1 = ((((((((-0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                          Wv8l2l1 + (0.25 * ((b_y0 - ye) * ((4.0 * x0 - 8.0 * xe) +
                                                            4.0) - ye * (4.0 * x0 - 4.0 * xe)) * (z0 - ze) - 0.25 * ze * (4.0 * x0 -
                                                                                                                          4.0 * xe) * (b_y0 - ye)) * y) + (0.33333333333333331 * ze * ((b_y0 - ye) *
                                                                                                                                                                                       ((4.0 * x0 - 8.0 * xe) + 4.0) - ye * (4.0 * x0 - 4.0 * xe))
                                                                                                                                                           + 0.33333333333333331 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * ((4.0 * x0
                                                                                                                                                                                                                            - 8.0 * xe) + 4.0)) * (z0 - ze)) * b_y) + (0.5 * ze * ((4.0 * xe - 4.0) *
                                                                                                                                                                                                                                                                                   (b_y0 - ye) + ye * ((4.0 * x0 - 8.0 * xe) + 4.0)) + 0.5 *
                                                                                                                                                                                                                                                                       ye * (4.0 * xe - 4.0) * (z0 - ze)) * (l1 * l1)) + ye * ze * (4.0 * xe -
                                                                                                                                                                                                                                                                                                                                    4.0) * l1) + 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) * c_y)
                     + (0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) - 0.25 *
                        ((b_y0 - ye) * ((4.0 * x0 - 8.0 * xe) + 4.0) - ye * (4.0 *
                                                                             x0 - 4.0 * xe)) * (z0 - ze)) * d_y) + (-0.33333333333333331 * ze * ((b_y0
                                                                                                                                                  - ye) * ((4.0 * x0 - 8.0 * xe) + 4.0) - ye * (4.0 * x0 - 4.0 * xe)) -
                                                                                                                    0.33333333333333331 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * ((4.0 * x0 -
                                                                                                                                                                                   8.0 * xe) + 4.0)) * (z0 - ze)) * e_y) + (-0.5 * ze * ((4.0 * xe - 4.0) *
                                                                                                                                                                                                                                         (b_y0 - ye) + ye * ((4.0 * x0 - 8.0 * xe) + 4.0)) - 0.5 * ye * (4.0 * xe -
                                                                                                                                                                                                                                                                                                         4.0) * (z0 - ze)) * (l2 * l2)) - ye * ze * (4.0 * xe - 4.0) * l2;

        /* 'find_weighting_trilin:31' Wv7l2l1 = - 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l1^5 + (1/4*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze))*(y0 - ye) - 1/4*ye*(4*z0 - 4*ze)*(x0 - xe))*l1^4 + (1/3*ye*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze)) + 1/3*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4))*(y0 - ye))*l1^3 + (1/2*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4)) + 1/2*xe*(4*ze - 4)*(y0 - ye))*l1^2 + xe*ye*(4*ze - 4)*l1 + 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l2^5 + (1/4*ye*(4*z0 - 4*ze)*(x0 - xe) - 1/4*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze))*(y0 - ye))*l2^4 + (- 1/3*ye*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze)) - 1/3*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4))*(y0 - ye))*l2^3 + (- 1/2*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4)) - 1/2*xe*(4*ze - 4)*(y0 - ye))*l2^2 - xe*ye*(4*ze - 4)*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv7l2l1 = ((((((((-0.2 * (4.0 * z0 - 4.0 * ze) * (x0 - xe) * (b_y0 - ye) *
                          Wv8l2l1 + (0.25 * ((x0 - xe) * ((4.0 * z0 - 8.0 * ze) +
                                                          4.0) - xe * (4.0 * z0 - 4.0 * ze)) * (b_y0 - ye) - 0.25 * ye * (4.0 * z0 -
                                                                                                                          4.0 * ze) * (x0 - xe)) * y) + (0.33333333333333331 * ye * ((x0 - xe) *
                                                                                                                                                                                     ((4.0 * z0 - 8.0 * ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze)) +
                                                                                                                                                         0.33333333333333331 * ((4.0 * ze - 4.0) * (x0 - xe) + xe * ((4.0 * z0 -
                                                                                                                                                                                                                      8.0 * ze) + 4.0)) * (b_y0 - ye)) * b_y) + (0.5 * ye * ((4.0 * ze - 4.0) *
                                                                                                                                                                                                                                                                             (x0 - xe) + xe * ((4.0 * z0 - 8.0 * ze) + 4.0)) + 0.5 * xe
                                                                                                                                                                                                                                                                 * (4.0 * ze - 4.0) * (b_y0 - ye)) * (l1 * l1)) + xe * ye * (4.0 * ze - 4.0)
                       * l1) + 0.2 * (4.0 * z0 - 4.0 * ze) * (x0 - xe) * (b_y0 - ye)
                      * c_y) + (0.25 * ye * (4.0 * z0 - 4.0 * ze) * (x0 - xe) - 0.25
                                * ((x0 - xe) * ((4.0 * z0 - 8.0 * ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze))
                                * (b_y0 - ye)) * d_y) + (-0.33333333333333331 * ye * ((x0 - xe) * ((4.0 *
                                                                                                    z0 - 8.0 * ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze)) - 0.33333333333333331 *
                                                         ((4.0 * ze - 4.0) * (x0 - xe) + xe * ((4.0 * z0 - 8.0 * ze) +
                                                                                               4.0)) * (b_y0 - ye)) * e_y) + (-0.5 * ye * ((4.0 * ze - 4.0) * (x0 - xe) +
                                                                                                                                           xe * ((4.0 * z0 - 8.0 * ze) + 4.0)) - 0.5 * xe * (4.0 * ze - 4.0) * (b_y0
                                                                                                                                                                                                                - ye)) * (l2 * l2)) - xe * ye * (4.0 * ze - 4.0) * l2;

        /* 'find_weighting_trilin:32' Wv8l2l1 = 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*ze*(4*x0 - 4*xe)*(y0 - ye) - 1/4*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe))*(z0 - ze))*l1^4 + (- 1/3*ze*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe)) - 1/3*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)))*l1^3 + (- 1/2*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)) - 2*xe*ye*(z0 - ze))*l1^2 - 4*xe*ye*ze*l1 - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe))*(z0 - ze) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l2^4 + (1/3*ze*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe)) + 1/3*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)))*l2^3 + (1/2*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)) + 2*xe*ye*(z0 - ze))*l2^2 + 4*xe*ye*ze*l2; */
        Wv8l2l1 = pow(l1, 5.0);
        y = pow(l1, 4.0);
        b_y = pow(l1, 3.0);
        c_y = pow(l2, 5.0);
        d_y = pow(l2, 4.0);
        e_y = pow(l2, 3.0);
        Wv8l2l1 = ((((((((0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                          Wv8l2l1 + (0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye)
                                     - 0.25 * ((4.0 * x0 - 8.0 * xe) * (b_y0 - ye) - ye * (4.0 * x0 - 4.0 * xe))
                                     * (z0 - ze)) * y) + (-0.33333333333333331 * ze * ((4.0 * x0 - 8.0 * xe) *
                                                                                       (b_y0 - ye) - ye * (4.0 * x0 - 4.0 * xe)) -
                                                          0.33333333333333331 * (z0 - ze) * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0
                                                                                                                            - 8.0 * xe))) * b_y) + (-0.5 * ze * (4.0 * xe * (b_y0 - ye) + ye * (4.0 *
                                                                                                                                                                                                x0 - 8.0 * xe)) - 2.0 * xe * ye * (z0 - ze)) * (l1 * l1)) - 4.0 * xe * ye *
                       ze * l1) - 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 -
                                                                               ze) * c_y) + (0.25 * ((4.0 * x0 - 8.0 * xe) * (b_y0 - ye) - ye * (4.0 * x0
                                                                                                                                                 - 4.0 * xe)) * (z0 - ze) - 0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye))
                     * d_y) + (0.33333333333333331 * ze * ((4.0 * x0 - 8.0 * xe) *
                                                           (b_y0 - ye) - ye * (4.0 * x0 - 4.0 * xe)) + 0.33333333333333331 * (z0 - ze)
                               * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0 - 8.0 * xe))) * e_y) + (0.5 *
                                                                                                  ze * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0 - 8.0 * xe)) + 2.0
                                                                                                  * xe * ye * (z0 - ze)) * (l2 * l2)) + 4.0 * xe * ye * ze * l2;
    } else {
        /* 'find_weighting_trilin:34' else */
        /*  l1 is less than 1/2, l2 is greater than 1/2 */
        /* 'find_weighting_trilin:36' tmp_l1 = l1; */
        /* 'find_weighting_trilin:37' tmp_l2 = l2; */
        /* 'find_weighting_trilin:39' l2 = 1/2; */
        /* 'find_weighting_trilin:40' l1 = tmp_l1; */
        /* 'find_weighting_trilin:42' Wv1l2l1 = 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) + (x0 - xe)*(y0 - ye)*(ze - 1))*l1^4 + (4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) + 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l1^3 + 2*(xe - 1)*(ye - 1)*(ze - 1)*l1^2 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- ((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - (x0 - xe)*(y0 - ye)*(ze - 1))*l2^4 + (- 4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) - 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l2^3 - 2*(xe - 1)*(ye - 1)*(ze - 1)*l2^2; */
        /* 'find_weighting_trilin:43' Wv2l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1))*(x0 - xe) - xe*(y0 - ye)*(z0 - ze))*l1^4 + (- 1/3*xe*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) - 4/3*(x0 - xe)*(ye - 1)*(ze - 1))*l1^3 - 2*xe*(ye - 1)*(ze - 1)*l1^2 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1))*(x0 - xe) + xe*(y0 - ye)*(z0 - ze))*l2^4 + (1/3*xe*(4*(y0 - ye)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) + 4/3*(x0 - xe)*(ye - 1)*(ze - 1))*l2^3 + 2*xe*(ye - 1)*(ze - 1)*l2^2; */
        /* 'find_weighting_trilin:44' Wv3l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1))*(y0 - ye) - ye*(x0 - xe)*(z0 - ze))*l1^4 + (- 1/3*ye*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) - 4/3*(y0 - ye)*(xe - 1)*(ze - 1))*l1^3 - 2*ye*(xe - 1)*(ze - 1)*l1^2 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1))*(y0 - ye) + ye*(x0 - xe)*(z0 - ze))*l2^4 + (1/3*ye*(4*(x0 - xe)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) + 4/3*(y0 - ye)*(xe - 1)*(ze - 1))*l2^3 + 2*ye*(xe - 1)*(ze - 1)*l2^2; */
        /* 'find_weighting_trilin:45' Wv4l2l1 = - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1))*(z0 - ze) - ze*(x0 - xe)*(y0 - ye))*l1^4 + (- 1/3*ze*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) - 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l1^3 - 2*ze*(xe - 1)*(ye - 1)*l1^2 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1))*(z0 - ze) + ze*(x0 - xe)*(y0 - ye))*l2^4 + (1/3*ze*(4*(x0 - xe)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) + 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l2^3 + 2*ze*(xe - 1)*(ye - 1)*l2^2; */
        /* 'find_weighting_trilin:46' Wv5l2l1 = 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l1^5 + (1/4*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye))*(z0 - ze) + 1/4*ze*(4*y0 - 4*ye)*(x0 - xe))*l1^4 + (1/3*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye)) + 1/3*xe*(4*ye - 4)*(z0 - ze))*l1^3 + 1/2*xe*ze*(4*ye - 4)*l1^2 - 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l2^5 + (- 1/4*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye))*(z0 - ze) - 1/4*ze*(4*y0 - 4*ye)*(x0 - xe))*l2^4 + (- 1/3*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 4*ye)) - 1/3*xe*(4*ye - 4)*(z0 - ze))*l2^3 - 1/2*xe*ze*(4*ye - 4)*l2^2; */
        /* 'find_weighting_trilin:47' Wv6l2l1 = 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe))*(z0 - ze) + 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l1^4 + (1/3*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe)) + 1/3*ye*(4*xe - 4)*(z0 - ze))*l1^3 + 1/2*ye*ze*(4*xe - 4)*l1^2 - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe))*(z0 - ze) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l2^4 + (- 1/3*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 4*xe)) - 1/3*ye*(4*xe - 4)*(z0 - ze))*l2^3 - 1/2*ye*ze*(4*xe - 4)*l2^2; */
        /* 'find_weighting_trilin:48' Wv7l2l1 = 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l1^5 + (1/4*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze))*(y0 - ye) + 1/4*ye*(4*z0 - 4*ze)*(x0 - xe))*l1^4 + (1/3*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze)) + 1/3*xe*(4*ze - 4)*(y0 - ye))*l1^3 + 1/2*xe*ye*(4*ze - 4)*l1^2 - 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l2^5 + (- 1/4*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze))*(y0 - ye) - 1/4*ye*(4*z0 - 4*ze)*(x0 - xe))*l2^4 + (- 1/3*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 4*ze)) - 1/3*xe*(4*ze - 4)*(y0 - ye))*l2^3 - 1/2*xe*ye*(4*ze - 4)*l2^2; */
        /* 'find_weighting_trilin:49' Wv8l2l1 = - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (- 1/4*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l1^4 + (- 1/3*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) - 4/3*xe*ye*(z0 - ze))*l1^3 - 2*xe*ye*ze*l1^2 + 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) + 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l2^4 + (1/3*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 4*xe)) + 4/3*xe*ye*(z0 - ze))*l2^3 + 2*xe*ye*ze*l2^2; */
        /* 'find_weighting_trilin:51' l1 = 1/2; */
        /* 'find_weighting_trilin:52' l2 = tmp_l2; */
        /* 'find_weighting_trilin:54' Wv1l2l1 = Wv1l2l1 + - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + ((x0 - xe)*(y0 - ye)*(z0 - ze) - ((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - (x0 - xe)*(y0 - ye)*(ze - 1))*l1^4 + (4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - 4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) + 4/3*(x0 - xe)*(y0 - ye)*(ze - 1) - 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l1^3 + (2*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) + 2*(z0 - ze)*(xe - 1)*(ye - 1) - 2*(xe - 1)*(ye - 1)*(ze - 1))*l1^2 + 4*(xe - 1)*(ye - 1)*(ze - 1)*l1 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - (x0 - xe)*(y0 - ye)*(z0 - ze) + (x0 - xe)*(y0 - ye)*(ze - 1))*l2^4 + (4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1) - 4/3*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(z0 - ze) - 4/3*(x0 - xe)*(y0 - ye)*(ze - 1) + 4/3*(z0 - ze)*(xe - 1)*(ye - 1))*l2^3 + (2*(xe - 1)*(ye - 1)*(ze - 1) - 2*(z0 - ze)*(xe - 1)*(ye - 1) - 2*((x0 - xe)*(ye - 1) + (y0 - ye)*(xe - 1))*(ze - 1))*l2^2 - 4*(xe - 1)*(ye - 1)*(ze - 1)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv1l2l1 = ((((((((((((((((0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y +
                                  (((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0))
                                   * (z0 - ze) + (x0 - xe) * (b_y0 - ye) * (ze - 1.0)) * d_y) +
                                 (1.3333333333333333 * ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) *
                                  (ze - 1.0) + 1.3333333333333333 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) *
                                 e_y) + 2.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0) * (l1 * l1)) - 0.8 * (x0
                                                                                                         - xe) * (b_y0 - ye) * (z0 - ze) * 0.03125) + (-((x0 - xe) * (ye - 1.0) +
                                                                                                                                                         (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) - (x0 - xe) * (b_y0 - ye) * (ze -
                                                                                                                                                                                                                            1.0)) * 0.0625) + (-1.3333333333333333 * ((x0 - xe) * (ye - 1.0) + (b_y0 -
                                                                                                                                                                                                                                                                                                ye) * (xe - 1.0)) * (ze - 1.0) - 1.3333333333333333 * (z0 - ze) * (xe -
                                                                                                                                                                                                                                                                                                                                                                   1.0) * (ye - 1.0)) * 0.125) - 2.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0) *
                            0.25) + -0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) *
                           0.03125) + (((x0 - xe) * (b_y0 - ye) * (z0 - ze) - ((x0 -
                                                                                xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0 - ze)) - (x0 - xe) *
                                       (b_y0 - ye) * (ze - 1.0)) * 0.0625) + (((1.3333333333333333 * ((x0 - xe) *
                                                                                                      (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) -
                                                                                1.3333333333333333 * ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) *
                                                                                (ze - 1.0)) + 1.3333333333333333 * (x0 - xe) * (b_y0 -
                                                                                                                                ye) * (ze - 1.0)) - 1.3333333333333333 * (z0 - ze) * (xe - 1.0) * (ye -
                                                                                                                                                                                                   1.0)) * 0.125) + ((2.0 * ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0))
                                                                                                                                                                                                                      * (ze - 1.0) + 2.0 * (z0 - ze) * (xe - 1.0) * (ye - 1.0))
                                                                                                                                                                                                                     - 2.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0)) * 0.25) +
                       4.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0) * 0.5) + 0.8 * (x0
                                                                                  - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1) + ((((x0 - xe) * (ye - 1.0) +
                                                                                                                                  (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) - (x0 - xe) * (b_y0 - ye) * (z0 - ze))
                                                                                                                                + (x0 - xe) * (b_y0 - ye) * (ze - 1.0)) * y) + (((1.3333333333333333 *
                                                                                                                                                                                  ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) * (ze - 1.0) -
                                                                                                                                                                                  1.3333333333333333 * ((x0 - xe) * (ye - 1.0) + (b_y0 - ye) * (xe - 1.0)) *
                                                                                                                                                                                  (z0 - ze)) - 1.3333333333333333 * (x0 - xe) * (b_y0 - ye) *
                                                                                                                                                                                 (ze - 1.0)) + 1.3333333333333333 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) *
                    b_y) + ((2.0 * (xe - 1.0) * (ye - 1.0) * (ze - 1.0) - 2.0 * (z0
                                                                                 - ze) * (xe - 1.0) * (ye - 1.0)) - 2.0 * ((x0 - xe) * (ye - 1.0) + (b_y0 -
                                                                                                                                                     ye) * (xe - 1.0)) * (ze - 1.0)) * (l2 * l2)) - 4.0 * (xe - 1.0) * (ye -
                                                                                                                                                                                                                        1.0) * (ze - 1.0) * l2;

        /* 'find_weighting_trilin:55' Wv2l2l1 = Wv2l2l1 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)) + xe*(y0 - ye)*(z0 - ze))*l1^4 + (1/3*xe*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)) - 1/3*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)))*l1^3 + (- 1/2*xe*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) - 2*(x0 - xe)*(ye - 1)*(ze - 1))*l1^2 - 4*xe*(ye - 1)*(ze - 1)*l1 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)) - xe*(y0 - ye)*(z0 - ze))*l2^4 + (1/3*(x0 - xe)*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) - 1/3*xe*(4*(y0 - ye)*(ze - 1) - 4*(y0 - ye)*(z0 - ze) + 4*(z0 - ze)*(ye - 1)))*l2^3 + (1/2*xe*(4*(y0 - ye)*(ze - 1) - 4*(ye - 1)*(ze - 1) + 4*(z0 - ze)*(ye - 1)) + 2*(x0 - xe)*(ye - 1)*(ze - 1))*l2^2 + 4*xe*(ye - 1)*(ze - 1)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv2l2l1 = ((((((((((((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y +
                                  (-0.25 * (4.0 * (b_y0 - ye) * (ze - 1.0) + 4.0 *
                                            (z0 - ze) * (ye - 1.0)) * (x0 - xe) - xe * (b_y0 - ye) * (z0 - ze)) * d_y)
                                 + (-0.33333333333333331 * xe * (4.0 * (b_y0 - ye) * (ze - 1.0) + 4.0 * (z0
                                                                                                         - ze) * (ye - 1.0)) - 1.3333333333333333 * (x0 - xe) * (ye - 1.0) * (ze -
                                                                                                                                                                              1.0)) * e_y) - 2.0 * xe * (ye - 1.0) * (ze - 1.0) * (l1 * l1)) + 0.8 * (x0
                                                                                                                                                                                                                                                      - xe) * (b_y0 - ye) * (z0 - ze) * 0.03125) + (0.25 * (4.0 * (b_y0 - ye) *
                                                                                                                                                                                                                                                                                                            (ze - 1.0) + 4.0 * (z0 - ze) * (ye - 1.0)) * (x0 -
                                                                                                                                                                                                                                                                                                                                                          xe) + xe * (b_y0 - ye) * (z0 - ze)) * 0.0625) + (0.33333333333333331 * xe *
                                                                                                                                                                                                                                                                                                                                                                                                           (4.0 * (b_y0 - ye) * (ze - 1.0) + 4.0 * (z0 - ze) *
                                                                                                                                                                                                                                                                                                                                                                                                            (ye - 1.0)) + 1.3333333333333333 * (x0 - xe) * (ye - 1.0) * (ze - 1.0)) *
                             0.125) + 2.0 * xe * (ye - 1.0) * (ze - 1.0) * 0.25) +
                           0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * 0.03125) +
                          (0.25 * (x0 - xe) * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0
                                                * (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0 - ze) * (ye - 1.0)) + xe * (b_y0 -
                                                                                                                   ye) * (z0 - ze)) * 0.0625) + (0.33333333333333331 * xe * ((4.0 * (b_y0 -
                                                                                                                                                                                     ye) * (ze - 1.0) - 4.0 * (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0 - ze) * (ye
                                                                                                                                                                                                                                                            - 1.0)) - 0.33333333333333331 * (x0 - xe) * ((4.0 * (b_y0 - ye) * (ze -
                                                                                                                                                                                                                                                                                                                               1.0) - 4.0 * (ye - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (ye - 1.0))) *
                         0.125) + (-0.5 * xe * ((4.0 * (b_y0 - ye) * (ze - 1.0) -
                                                 4.0 * (ye - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (ye - 1.0)) - 2.0 * (x0
                                                                                                                         - xe) * (ye - 1.0) * (ze - 1.0)) * 0.25) - 4.0 * xe * (ye - 1.0) * (ze -
                                                                                                                                                                                             1.0) * 0.5) - 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * Wv8l2l1) +
                     (-0.25 * (x0 - xe) * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0 *
                                            (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0 - ze) * (ye - 1.0)) - xe * (b_y0 - ye)
                      * (z0 - ze)) * y) + (0.33333333333333331 * (x0 - xe) * ((4.0 *
                                                                               (b_y0 - ye) * (ze - 1.0) - 4.0 * (ye - 1.0) * (ze - 1.0)) +
                                                                              4.0 * (z0 - ze) * (ye - 1.0)) - 0.33333333333333331 * xe * ((4.0 * (b_y0 -
                                                                                                                                                  ye) * (ze - 1.0) - 4.0 * (b_y0 - ye) * (z0 - ze)) + 4.0 * (z0 - ze) * (ye
                                                                                                                                                                                                                         - 1.0))) * b_y) + (0.5 * xe * ((4.0 * (b_y0 - ye) * (ze - 1.0) - 4.0 * (ye
                                                                                                                                                                                                                                                                                                 - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (ye - 1.0)) + 2.0 * (x0 - xe) *
                                                                                                                                                                                                                                            (ye - 1.0) * (ze - 1.0)) * (l2 * l2)) + 4.0 * xe * (ye
                                                                                                                                                                                                                                                                                                - 1.0) * (ze - 1.0) * l2;

        /* 'find_weighting_trilin:56' Wv3l2l1 = Wv3l2l1 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)) + ye*(x0 - xe)*(z0 - ze))*l1^4 + (1/3*ye*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)) - 1/3*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)))*l1^3 + (- 1/2*ye*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) - 2*(y0 - ye)*(xe - 1)*(ze - 1))*l1^2 - 4*ye*(xe - 1)*(ze - 1)*l1 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)) - ye*(x0 - xe)*(z0 - ze))*l2^4 + (1/3*(y0 - ye)*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) - 1/3*ye*(4*(x0 - xe)*(ze - 1) - 4*(x0 - xe)*(z0 - ze) + 4*(z0 - ze)*(xe - 1)))*l2^3 + (1/2*ye*(4*(x0 - xe)*(ze - 1) - 4*(xe - 1)*(ze - 1) + 4*(z0 - ze)*(xe - 1)) + 2*(y0 - ye)*(xe - 1)*(ze - 1))*l2^2 + 4*ye*(xe - 1)*(ze - 1)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv3l2l1 = ((((((((((((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y +
                                  (-0.25 * (4.0 * (x0 - xe) * (ze - 1.0) + 4.0 * (z0
                                                                                  - ze) * (xe - 1.0)) * (b_y0 - ye) - ye * (x0 - xe) * (z0 - ze)) * d_y) + (
                                     -0.33333333333333331 * ye * (4.0 * (x0 - xe) * (ze - 1.0) + 4.0 * (z0 - ze)
                                                                  * (xe - 1.0)) - 1.3333333333333333 * (b_y0 - ye) * (xe - 1.0) * (ze - 1.0))
                                 * e_y) - 2.0 * ye * (xe - 1.0) * (ze - 1.0) * (l1 * l1)) + 0.8 * (x0 - xe)
                               * (b_y0 - ye) * (z0 - ze) * 0.03125) + (0.25 * (4.0 * (x0 - xe) * (ze -
                                                                                                  1.0) + 4.0 * (z0 - ze) * (xe - 1.0)) * (b_y0 - ye) + ye * (x0 - xe) * (z0
                                                                                                                                                                         - ze)) * 0.0625) + (0.33333333333333331 * ye * (4.0 * (x0 - xe) * (ze -
                                                                                                                                                                                                                                            1.0) + 4.0 * (z0 - ze) * (xe - 1.0)) + 1.3333333333333333 * (b_y0 - ye) *
                                                                                                                                                                                             (xe - 1.0) * (ze - 1.0)) * 0.125) + 2.0 * ye * (xe -
                                                                                                                                                                                                                                             1.0) * (ze - 1.0) * 0.25) + 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) *
                           0.03125) + (0.25 * (b_y0 - ye) * ((4.0 * (x0 - xe) * (ze
                                                                                 - 1.0) - 4.0 * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) * (xe - 1.0)) + ye
                                       * (x0 - xe) * (z0 - ze)) * 0.0625) + (0.33333333333333331 * ye * ((4.0 *
                                                                                                          (x0 - xe) * (ze - 1.0) - 4.0 * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) *
                                                                                                         (xe - 1.0)) - 0.33333333333333331 * (b_y0 - ye) * ((4.0 *
                                                                                                                                                             (x0 - xe) * (ze - 1.0) - 4.0 * (xe - 1.0) * (ze - 1.0))
                                                                                                                                                            + 4.0 * (z0 - ze) * (xe - 1.0))) * 0.125) + (-0.5 * ye * ((4.0 * (x0 - xe)
                                                                                                                                                                                                                       * (ze - 1.0) - 4.0 * (xe - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (xe -
                                                                                                                                                                                                                                                                                          1.0)) - 2.0 * (b_y0 - ye) * (xe - 1.0) * (ze - 1.0)) * 0.25) - 4.0 * ye *
                       (xe - 1.0) * (ze - 1.0) * 0.5) - 0.8 * (x0 - xe) * (b_y0 - ye)
                      * (z0 - ze) * Wv8l2l1) + (-0.25 * (b_y0 - ye) * ((4.0 * (x0 -
                                                                               xe) * (ze - 1.0) - 4.0 * (x0 - xe) * (z0 - ze)) + 4.0 * (z0 - ze) * (xe -
                                                                                                                                                    1.0)) - ye * (x0 - xe) * (z0 - ze)) * y) + (0.33333333333333331 * (b_y0 -
                                                                                                                                                                                                                       ye) * ((4.0 * (x0 - xe) * (ze - 1.0) - 4.0 * (xe - 1.0) * (ze - 1.0)) +
                                                                                                                                                                                                                              4.0 * (z0 - ze) * (xe - 1.0)) - 0.33333333333333331 * ye * ((4.0 *
                                                                                                                                                                                                                                                                                           (x0 - xe) * (ze - 1.0) - 4.0 * (x0 - xe) * (z0 - ze)) + 4.0 *
                                                                                                                                                                                                                                                                                          (z0 - ze) * (xe - 1.0))) * b_y) + (0.5 * ye * ((4.0 * (x0 - xe)
                                                                                                                                                                                                                                                                                                                                          * (ze - 1.0) - 4.0 * (xe - 1.0) * (ze - 1.0)) + 4.0 * (z0 - ze) * (xe -
                                                                                                                                                                                                                                                                                                                                                                                                             1.0)) + 2.0 * (b_y0 - ye) * (xe - 1.0) * (ze - 1.0)) * (l2 * l2)) + 4.0 *
                ye * (xe - 1.0) * (ze - 1.0) * l2;

        /* 'find_weighting_trilin:57' Wv4l2l1 = Wv4l2l1 + 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)) + ze*(x0 - xe)*(y0 - ye))*l1^4 + (1/3*ze*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)) - 1/3*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)))*l1^3 + (- 1/2*ze*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) - 2*(z0 - ze)*(xe - 1)*(ye - 1))*l1^2 - 4*ze*(xe - 1)*(ye - 1)*l1 - 4/5*(x0 - xe)*(y0 - ye)*(z0 - ze)*l2^5 + (- 1/4*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)) - ze*(x0 - xe)*(y0 - ye))*l2^4 + (1/3*(z0 - ze)*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) - 1/3*ze*(4*(x0 - xe)*(ye - 1) - 4*(x0 - xe)*(y0 - ye) + 4*(y0 - ye)*(xe - 1)))*l2^3 + (1/2*ze*(4*(x0 - xe)*(ye - 1) - 4*(xe - 1)*(ye - 1) + 4*(y0 - ye)*(xe - 1)) + 2*(z0 - ze)*(xe - 1)*(ye - 1))*l2^2 + 4*ze*(xe - 1)*(ye - 1)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv4l2l1 = ((((((((((((((((-0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) * c_y +
                                  (-0.25 * (4.0 * (x0 - xe) * (ye - 1.0) + 4.0 *
                                            (b_y0 - ye) * (xe - 1.0)) * (z0 - ze) - ze * (x0 - xe) * (b_y0 - ye)) *
                                  d_y) + (-0.33333333333333331 * ze * (4.0 * (x0 - xe) * (ye - 1.0) + 4.0 *
                                                                       (b_y0 - ye) * (xe - 1.0)) - 1.3333333333333333 *
                                          (z0 - ze) * (xe - 1.0) * (ye - 1.0)) * e_y) - 2.0 * ze * (xe - 1.0)
                                * (ye - 1.0) * (l1 * l1)) + 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) *
                               0.03125) + (0.25 * (4.0 * (x0 - xe) * (ye - 1.0) + 4.0 * (b_y0 - ye) * (xe
                                                                                                       - 1.0)) * (z0 - ze) + ze * (x0 - xe) * (b_y0 - ye)) * 0.0625) +
                             (0.33333333333333331 * ze * (4.0 * (x0 - xe) * (ye -
                                                                             1.0) + 4.0 * (b_y0 - ye) * (xe - 1.0)) + 1.3333333333333333 * (z0 - ze) *
                              (xe - 1.0) * (ye - 1.0)) * 0.125) + 2.0 * ze * (xe -
                                                                              1.0) * (ye - 1.0) * 0.25) + 0.8 * (x0 - xe) * (b_y0 - ye) * (z0 - ze) *
                           0.03125) + (0.25 * (z0 - ze) * ((4.0 * (x0 - xe) * (ye -
                                                                               1.0) - 4.0 * (x0 - xe) * (b_y0 - ye)) + 4.0 * (b_y0 - ye) * (xe - 1.0)) +
                                       ze * (x0 - xe) * (b_y0 - ye)) * 0.0625) + (0.33333333333333331 * ze *
                                                                                  ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (x0 - xe) * (b_y0 - ye)) + 4.0 *
                                                                                   (b_y0 - ye) * (xe - 1.0)) - 0.33333333333333331 * (z0 - ze) * ((4.0 * (x0
                                                                                                                                                          - xe) * (ye - 1.0) - 4.0 * (xe - 1.0) * (ye - 1.0)) + 4.0 * (b_y0 - ye) *
                                                                                                                                                  (xe - 1.0))) * 0.125) + (-0.5 * ze * ((4.0 * (x0 - xe) *
                                                                                                                                                                                         (ye - 1.0) - 4.0 * (xe - 1.0) * (ye - 1.0)) + 4.0 * (b_y0
                                                                                                                                                                                                                                              - ye) * (xe - 1.0)) - 2.0 * (z0 - ze) * (xe - 1.0) * (ye - 1.0)) * 0.25) -
                       4.0 * ze * (xe - 1.0) * (ye - 1.0) * 0.5) - 0.8 * (x0 - xe) *
                      (b_y0 - ye) * (z0 - ze) * Wv8l2l1) + (-0.25 * (z0 - ze) *
                                                            ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (x0 - xe) * (b_y0 - ye)) + 4.0 *
                                                             (b_y0 - ye) * (xe - 1.0)) - ze * (x0 - xe) * (b_y0 - ye)) * y) +
                    (0.33333333333333331 * (z0 - ze) * ((4.0 * (x0 - xe) * (ye - 1.0)
                                                         - 4.0 * (xe - 1.0) * (ye - 1.0)) + 4.0 * (b_y0 - ye) * (xe - 1.0)) -
                     0.33333333333333331 * ze * ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0
                                                  * (x0 - xe) * (b_y0 - ye)) + 4.0 * (b_y0 - ye) * (xe - 1.0))) * b_y) +
                   (0.5 * ze * ((4.0 * (x0 - xe) * (ye - 1.0) - 4.0 * (xe - 1.0) *
                                 (ye - 1.0)) + 4.0 * (b_y0 - ye) * (xe - 1.0)) + 2.0 * (z0 - ze) * (xe -
                                                                                                    1.0) * (ye - 1.0)) * (l2 * l2)) + 4.0 * ze * (xe - 1.0) * (ye - 1.0) * l2;

        /* 'find_weighting_trilin:58' Wv5l2l1 = Wv5l2l1 - 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l1^5 + (1/4*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye))*(z0 - ze) - 1/4*ze*(4*y0 - 4*ye)*(x0 - xe))*l1^4 + (1/3*ze*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye)) + 1/3*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4))*(z0 - ze))*l1^3 + (1/2*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4)) + 1/2*xe*(4*ye - 4)*(z0 - ze))*l1^2 + xe*ze*(4*ye - 4)*l1 + 1/5*(4*y0 - 4*ye)*(x0 - xe)*(z0 - ze)*l2^5 + (1/4*ze*(4*y0 - 4*ye)*(x0 - xe) - 1/4*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye))*(z0 - ze))*l2^4 + (- 1/3*ze*((x0 - xe)*(4*y0 - 8*ye + 4) - xe*(4*y0 - 4*ye)) - 1/3*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4))*(z0 - ze))*l2^3 + (- 1/2*ze*((4*ye - 4)*(x0 - xe) + xe*(4*y0 - 8*ye + 4)) - 1/2*xe*(4*ye - 4)*(z0 - ze))*l2^2 - xe*ze*(4*ye - 4)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv5l2l1 = ((((((((((((((((0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 -
                                                                               ze) * c_y + (0.25 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * (4.0 * b_y0 - 4.0
                                                                                                                                         * ye)) * (z0 - ze) + 0.25 * ze * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe)) *
                                  d_y) + (0.33333333333333331 * ze * ((4.0 * ye - 4.0) * (x0 - xe) + xe *
                                                                      (4.0 * b_y0 - 4.0 * ye)) + 0.33333333333333331 * xe * (4.0 * ye - 4.0) *
                                          (z0 - ze)) * e_y) + 0.5 * xe * ze * (4.0 * ye - 4.0) * (l1 * l1))
                               - 0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 - ze) * 0.03125) +
                              (-0.25 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * (4.0 *
                                                                             b_y0 - 4.0 * ye)) * (z0 - ze) - 0.25 * ze * (4.0 * b_y0 - 4.0 * ye) * (x0
                                                                                                                                                    - xe)) * 0.0625) + (-0.33333333333333331 * ze * ((4.0 * ye - 4.0) * (x0 -
                                                                                                                                                                                                                         xe) + xe * (4.0 * b_y0 - 4.0 * ye)) - 0.33333333333333331 * xe * (4.0 * ye
                                                                                                                                                                                                                                                                                           - 4.0) * (z0 - ze)) * 0.125) - 0.5 * xe * ze * (4.0 * ye - 4.0) * 0.25) -
                           0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 - ze) *
                           0.03125) + (0.25 * ((x0 - xe) * ((4.0 * b_y0 - 8.0 * ye)
                                                            + 4.0) - xe * (4.0 * b_y0 - 4.0 * ye)) * (z0 - ze) - 0.25 * ze * (4.0 *
                                                                                                                              b_y0 - 4.0 * ye) * (x0 - xe)) * 0.0625) + (0.33333333333333331 * ze * ((x0
                                                                                                                                                                                                      - xe) * ((4.0 * b_y0 - 8.0 * ye) + 4.0) - xe * (4.0 * b_y0 - 4.0 * ye)) +
                                                                                                                                                                         0.33333333333333331 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * ((4.0 * b_y0 -
                                                                                                                                                                                                                                      8.0 * ye) + 4.0)) * (z0 - ze)) * 0.125) + (0.5 * ze * ((4.0 * ye - 4.0) *
                                                                                                                                                                                                                                                                                             (x0 - xe) + xe * ((4.0 * b_y0 - 8.0 * ye) + 4.0)) + 0.5 *
                                                                                                                                                                                                                                                                                 xe * (4.0 * ye - 4.0) * (z0 - ze)) * 0.25) + xe * ze * (4.0 * ye - 4.0) *
                       0.5) + 0.2 * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) * (z0 - ze) *
                      Wv8l2l1) + (0.25 * ze * (4.0 * b_y0 - 4.0 * ye) * (x0 - xe) -
                                  0.25 * ((x0 - xe) * ((4.0 * b_y0 - 8.0 * ye) + 4.0) - xe * (4.0 * b_y0 -
                                                                                              4.0 * ye)) * (z0 - ze)) * y) + (-0.33333333333333331 * ze * ((x0 - xe) *
                                                                                                                                                           ((4.0 * b_y0 - 8.0 * ye) + 4.0) - xe * (4.0 * b_y0 - 4.0 * ye))
                                                                                                                              - 0.33333333333333331 * ((4.0 * ye - 4.0) * (x0 - xe) + xe * ((4.0 * b_y0
                                                                                                                                                                                             - 8.0 * ye) + 4.0)) * (z0 - ze)) * b_y) + (-0.5 * ze * ((4.0 * ye - 4.0) *
                                                                                                                                                                                                                                                     (x0 - xe) + xe * ((4.0 * b_y0 - 8.0 * ye) + 4.0)) - 0.5 * xe *
                                                                                                                                                                                                                                        (4.0 * ye - 4.0) * (z0 - ze)) * (l2 * l2)) - xe * ze * (4.0 * ye
                                                                                                                                                                                                                                                                                                - 4.0) * l2;

        /* 'find_weighting_trilin:59' Wv6l2l1 = Wv6l2l1 - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe))*(z0 - ze) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l1^4 + (1/3*ze*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe)) + 1/3*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4))*(z0 - ze))*l1^3 + (1/2*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4)) + 1/2*ye*(4*xe - 4)*(z0 - ze))*l1^2 + ye*ze*(4*xe - 4)*l1 + 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*ze*(4*x0 - 4*xe)*(y0 - ye) - 1/4*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe))*(z0 - ze))*l2^4 + (- 1/3*ze*((y0 - ye)*(4*x0 - 8*xe + 4) - ye*(4*x0 - 4*xe)) - 1/3*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4))*(z0 - ze))*l2^3 + (- 1/2*ze*((4*xe - 4)*(y0 - ye) + ye*(4*x0 - 8*xe + 4)) - 1/2*ye*(4*xe - 4)*(z0 - ze))*l2^2 - ye*ze*(4*xe - 4)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv6l2l1 = ((((((((((((((((0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 -
                                                                               ze) * c_y + (0.25 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * (4.0 * x0 - 4.0
                                                                                                                                           * xe)) * (z0 - ze) + 0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye)) *
                                  d_y) + (0.33333333333333331 * ze * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye *
                                                                      (4.0 * x0 - 4.0 * xe)) + 0.33333333333333331 * ye
                                          * (4.0 * xe - 4.0) * (z0 - ze)) * e_y) + 0.5 * ye * ze * (4.0 * xe
                                                                                                    - 4.0) * (l1 * l1)) - 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze)
                               * 0.03125) + (-0.25 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * (4.0 * x0 -
                                                                                             4.0 * xe)) * (z0 - ze) - 0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye)) *
                              0.0625) + (-0.33333333333333331 * ze * ((4.0 * xe -
                                                                       4.0) * (b_y0 - ye) + ye * (4.0 * x0 - 4.0 * xe)) - 0.33333333333333331 *
                                         ye * (4.0 * xe - 4.0) * (z0 - ze)) * 0.125) - 0.5 * ye * ze * (4.0 * xe -
                                                                                                        4.0) * 0.25) - 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                           0.03125) + (0.25 * ((b_y0 - ye) * ((4.0 * x0 - 8.0 * xe)
                                                              + 4.0) - ye * (4.0 * x0 - 4.0 * xe)) * (z0 - ze) - 0.25 * ze * (4.0 * x0 -
                                                                                                                              4.0 * xe) * (b_y0 - ye)) * 0.0625) + (0.33333333333333331 * ze * ((b_y0 -
                                                                                                                                                                                                 ye) * ((4.0 * x0 - 8.0 * xe) + 4.0) - ye * (4.0 * x0 - 4.0 * xe)) +
                                                                                                                                                                    0.33333333333333331 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * ((4.0 * x0 -
                                                                                                                                                                                                                                   8.0 * xe) + 4.0)) * (z0 - ze)) * 0.125) + (0.5 * ze * ((4.0 * xe - 4.0) *
                                                                                                                                                                                                                                                                                          (b_y0 - ye) + ye * ((4.0 * x0 - 8.0 * xe) + 4.0)) + 0.5 *
                                                                                                                                                                                                                                                                              ye * (4.0 * xe - 4.0) * (z0 - ze)) * 0.25) + ye * ze * (4.0 * xe - 4.0) *
                       0.5) + 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                      Wv8l2l1) + (0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) -
                                  0.25 * ((b_y0 - ye) * ((4.0 * x0 - 8.0 * xe) + 4.0) - ye * (4.0 * x0 - 4.0
                                                                                              * xe)) * (z0 - ze)) * y) + (-0.33333333333333331 * ze * ((b_y0 - ye) *
                                                                                                                                                       ((4.0 * x0 - 8.0 * xe) + 4.0) - ye * (4.0 * x0 - 4.0 * xe)) -
                                                                                                                          0.33333333333333331 * ((4.0 * xe - 4.0) * (b_y0 - ye) + ye * ((4.0 * x0 -
                                                                                                                                                                                         8.0 * xe) + 4.0)) * (z0 - ze)) * b_y) + (-0.5 * ze * ((4.0 * xe - 4.0) *
                                                                                                                                                                                                                                               (b_y0 - ye) + ye * ((4.0 * x0 - 8.0 * xe) + 4.0)) - 0.5 * ye * (4.0 * xe -
                                                                                                                                                                                                                                                                                                               4.0) * (z0 - ze)) * (l2 * l2)) - ye * ze * (4.0 * xe - 4.0) * l2;

        /* 'find_weighting_trilin:60' Wv7l2l1 = Wv7l2l1 - 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l1^5 + (1/4*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze))*(y0 - ye) - 1/4*ye*(4*z0 - 4*ze)*(x0 - xe))*l1^4 + (1/3*ye*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze)) + 1/3*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4))*(y0 - ye))*l1^3 + (1/2*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4)) + 1/2*xe*(4*ze - 4)*(y0 - ye))*l1^2 + xe*ye*(4*ze - 4)*l1 + 1/5*(4*z0 - 4*ze)*(x0 - xe)*(y0 - ye)*l2^5 + (1/4*ye*(4*z0 - 4*ze)*(x0 - xe) - 1/4*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze))*(y0 - ye))*l2^4 + (- 1/3*ye*((x0 - xe)*(4*z0 - 8*ze + 4) - xe*(4*z0 - 4*ze)) - 1/3*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4))*(y0 - ye))*l2^3 + (- 1/2*ye*((4*ze - 4)*(x0 - xe) + xe*(4*z0 - 8*ze + 4)) - 1/2*xe*(4*ze - 4)*(y0 - ye))*l2^2 - xe*ye*(4*ze - 4)*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv7l2l1 = ((((((((((((((((0.2 * (4.0 * z0 - 4.0 * ze) * (x0 - xe) * (b_y0 -
                                                                             ye) * c_y + (0.25 * ((4.0 * ze - 4.0) * (x0 - xe) + xe * (4.0 * z0 - 4.0 *
                                                                                                                                       ze)) * (b_y0 - ye) + 0.25 * ye * (4.0 * z0 - 4.0 * ze) * (x0 - xe)) * d_y)
                                 + (0.33333333333333331 * ye * ((4.0 * ze - 4.0) * (x0 - xe) + xe * (4.0 *
                                                                                                     z0 - 4.0 * ze)) + 0.33333333333333331 * xe * (4.0 * ze - 4.0) * (b_y0 - ye))
                                 * e_y) + 0.5 * xe * ye * (4.0 * ze - 4.0) * (l1 * l1)) - 0.2 * (4.0 * z0 -
                                                                                                 4.0 * ze) * (x0 - xe) * (b_y0 - ye) * 0.03125) + (-0.25 * ((4.0 * ze - 4.0)
                                                                                                                                                            * (x0 - xe) + xe * (4.0 * z0 - 4.0 * ze)) * (b_y0 - ye) - 0.25 * ye * (4.0
                                                                                                                                                                                                                                   * z0 - 4.0 * ze) * (x0 - xe)) * 0.0625) + (-0.33333333333333331 * ye *
                                                                                                                                                                                                                                                                              ((4.0 * ze - 4.0) * (x0 - xe) + xe * (4.0 * z0 - 4.0 * ze)) -
                                                                                                                                                                                                                                                                              0.33333333333333331 * xe * (4.0 * ze - 4.0) * (b_y0 - ye)) * 0.125) - 0.5 *
                            xe * ye * (4.0 * ze - 4.0) * 0.25) - 0.2 * (4.0 * z0 -
                                                                        4.0 * ze) * (x0 - xe) * (b_y0 - ye) * 0.03125) + (0.25 * ((x0 - xe) *
                                                                                                                                  ((4.0 * z0 - 8.0 * ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze)) * (b_y0 - ye)
                                                                                                                          - 0.25 * ye * (4.0 * z0 - 4.0 * ze) * (x0 - xe)) * 0.0625) +
                         (0.33333333333333331 * ye * ((x0 - xe) * ((4.0 * z0 - 8.0 *
                                                                    ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze)) + 0.33333333333333331 * ((4.0 *
                                                                                                                                       ze - 4.0) * (x0 - xe) + xe * ((4.0 * z0 - 8.0 * ze) + 4.0)) * (b_y0 - ye))
                         * 0.125) + (0.5 * ye * ((4.0 * ze - 4.0) * (x0 - xe) + xe *
                                                 ((4.0 * z0 - 8.0 * ze) + 4.0)) + 0.5 * xe * (4.0 * ze -
                                                                                              4.0) * (b_y0 - ye)) * 0.25) + xe * ye * (4.0 * ze - 4.0) * 0.5) + 0.2 *
                      (4.0 * z0 - 4.0 * ze) * (x0 - xe) * (b_y0 - ye) * Wv8l2l1) +
                     (0.25 * ye * (4.0 * z0 - 4.0 * ze) * (x0 - xe) - 0.25 * ((x0 -
                                                                               xe) * ((4.0 * z0 - 8.0 * ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze)) * (b_y0
                                                                                                                                                    - ye)) * y) + (-0.33333333333333331 * ye * ((x0 - xe) * ((4.0 * z0 - 8.0 *
                                                                                                                                                                                                              ze) + 4.0) - xe * (4.0 * z0 - 4.0 * ze)) - 0.33333333333333331 * ((4.0 *
                                                                                                                                                                                                                                                                                 ze - 4.0) * (x0 - xe) + xe * ((4.0 * z0 - 8.0 * ze) + 4.0)) * (b_y0 - ye))
                    * b_y) + (-0.5 * ye * ((4.0 * ze - 4.0) * (x0 - xe) + xe * ((4.0
                                                                                 * z0 - 8.0 * ze) + 4.0)) - 0.5 * xe * (4.0 * ze - 4.0) * (b_y0 - ye)) *
                   (l2 * l2)) - xe * ye * (4.0 * ze - 4.0) * l2;

        /* 'find_weighting_trilin:61' Wv8l2l1 = Wv8l2l1 + 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l1^5 + (1/4*ze*(4*x0 - 4*xe)*(y0 - ye) - 1/4*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe))*(z0 - ze))*l1^4 + (- 1/3*ze*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe)) - 1/3*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)))*l1^3 + (- 1/2*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)) - 2*xe*ye*(z0 - ze))*l1^2 - 4*xe*ye*ze*l1 - 1/5*(4*x0 - 4*xe)*(y0 - ye)*(z0 - ze)*l2^5 + (1/4*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe))*(z0 - ze) - 1/4*ze*(4*x0 - 4*xe)*(y0 - ye))*l2^4 + (1/3*ze*((4*x0 - 8*xe)*(y0 - ye) - ye*(4*x0 - 4*xe)) + 1/3*(z0 - ze)*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)))*l2^3 + (1/2*ze*(4*xe*(y0 - ye) + ye*(4*x0 - 8*xe)) + 2*xe*ye*(z0 - ze))*l2^2 + 4*xe*ye*ze*l2; */
        Wv8l2l1 = pow(l2, 5.0);
        y = pow(l2, 4.0);
        b_y = pow(l2, 3.0);
        c_y = pow(l1, 5.0);
        d_y = pow(l1, 4.0);
        e_y = pow(l1, 3.0);
        Wv8l2l1 = ((((((((((((((((-0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 -
                                                                                ze) * c_y + (-0.25 * (z0 - ze) * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0
                                                                                                                                                 - 4.0 * xe)) - 0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye)) * d_y) + (
                                     -0.33333333333333331 * ze * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0 - 4.0
                                                                                                 * xe)) - 1.3333333333333333 * xe * ye * (z0 - ze)) * e_y) - 2.0 * xe * ye *
                                ze * (l1 * l1)) + 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                               0.03125) + (0.25 * (z0 - ze) * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0 -
                                                                                              4.0 * xe)) + 0.25 * ze * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye)) * 0.0625) +
                             (0.33333333333333331 * ze * (4.0 * xe * (b_y0 - ye) +
                                                          ye * (4.0 * x0 - 4.0 * xe)) + 1.3333333333333333 * xe * ye * (z0 - ze)) *
                             0.125) + 2.0 * xe * ye * ze * 0.25) + 0.2 * (4.0 * x0 -
                                                                          4.0 * xe) * (b_y0 - ye) * (z0 - ze) * 0.03125) + (0.25 * ze * (4.0 * x0 -
                                                                                                                                         4.0 * xe) * (b_y0 - ye) - 0.25 * ((4.0 * x0 - 8.0 * xe) * (b_y0 - ye) - ye
                                                                                                                                                                           * (4.0 * x0 - 4.0 * xe)) * (z0 - ze)) * 0.0625) + (-0.33333333333333331 *
                                                                                                                                                                                                                              ze * ((4.0 * x0 - 8.0 * xe) * (b_y0 - ye) - ye * (4.0 * x0 - 4.0 * xe)) -
                                                                                                                                                                                                                              0.33333333333333331 * (z0 - ze) * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0
                                                                                                                                                                                                                                                                                                - 8.0 * xe))) * 0.125) + (-0.5 * ze * (4.0 * xe * (b_y0 - ye) + ye * (4.0 *
                                                                                                                                                                                                                                                                                                                                                                      x0 - 8.0 * xe)) - 2.0 * xe * ye * (z0 - ze)) * 0.25) - 4.0 * xe * ye * ze *
                       0.5) - 0.2 * (4.0 * x0 - 4.0 * xe) * (b_y0 - ye) * (z0 - ze) *
                      Wv8l2l1) + (0.25 * ((4.0 * x0 - 8.0 * xe) * (b_y0 - ye) - ye *
                                          (4.0 * x0 - 4.0 * xe)) * (z0 - ze) - 0.25 * ze * (4.0 * x0 -
                                                                                            4.0 * xe) * (b_y0 - ye)) * y) + (0.33333333333333331 * ze * ((4.0 * x0 -
                                                                                                                                                          8.0 * xe) * (b_y0 - ye) - ye * (4.0 * x0 - 4.0 * xe)) +
                                                                                                                             0.33333333333333331 * (z0 - ze) * (4.0 * xe * (b_y0 - ye) + ye * (4.0 * x0
                                                                                                                                                                                               - 8.0 * xe))) * b_y) + (0.5 * ze * (4.0 * xe * (b_y0 - ye) + ye * (4.0 *
                                                                                                                                                                                                                                                                  x0 - 8.0 * xe)) + 2.0 * xe * ye * (z0 - ze)) * (l2 * l2)) + 4.0 * xe * ye *
                ze * l2;
    }

    /*  old find weighting function */
    /*  if x1 > x2 */
    /*     error('x1 > x2')  */
    /*  end */
    /*   */
    /*  if x1 < -(depth_cell_size/2 + overhang) - eps(x1) */
    /*     x1 = -(depth_cell_size/2 + overhang); */
    /*  %     warning('x1 too low x1max = %0.3f x1 = %0.3f',-(depth_cell_size/2 + overhang) - eps(x1),x1)  */
    /*  elseif x2 > (depth_cell_size/2 + overhang) + eps(x2) */
    /*      x2 = (depth_cell_size/2 + overhang); */
    /*  %     warning('x2 too high x2max = %0.3f x2 = %0.3f',(depth_cell_size/2 + overhang) + eps(x2),x2 ) */
    /*       */
    /*  end */
    /*   */
    /*  if (x1 <= 0) && (x2 >= 0) */
    /*      A1 = 0.5 * (overhang + depth_cell_size/2 + x1)^2/(overhang + depth_cell_size/2)^2; */
    /*      A3 = 0.5 * (overhang + depth_cell_size/2 - x2)^2/(overhang + depth_cell_size/2)^2; */
    /*      w = 1 - A1 - A3; */
    /*  elseif (x1 <= 0) && (x2 <= 0) */
    /*      A1 = 0.5 * (overhang + depth_cell_size/2 + x1)^2/(overhang + depth_cell_size/2)^2; */
    /*      A2 = 0.5 * (overhang + depth_cell_size/2 + x2)^2/(overhang + depth_cell_size/2)^2; */
    /*      w = A2 - A1; */
    /*  elseif (x1 >= 0) && (x2 >= 0) */
    /*      A2 = 0.5 * (overhang + depth_cell_size/2 - x1)^2/(overhang + depth_cell_size/2)^2; */
    /*      A3 = 0.5 * (overhang + depth_cell_size/2 - x2)^2/(overhang + depth_cell_size/2)^2; */
    /*      w = A2 - A3; */
    /*  end */
    /*  if w < 0 */
    /*     error('w is negative') */
    /*  end */
    /* 'find_weighting_trilin:100' w = [Wv1l2l1 ; */
    /* 'find_weighting_trilin:101'     Wv2l2l1 ; */
    /* 'find_weighting_trilin:102'     Wv3l2l1 ; */
    /* 'find_weighting_trilin:103'     Wv4l2l1 ; */
    /* 'find_weighting_trilin:104'     Wv5l2l1 ; */
    /* 'find_weighting_trilin:105'     Wv6l2l1 ; */
    /* 'find_weighting_trilin:106'     Wv7l2l1 ; */
    /* 'find_weighting_trilin:107'     Wv8l2l1]; */
    w[0] = Wv1l2l1;
    w[1] = Wv2l2l1;
    w[2] = Wv3l2l1;
    w[3] = Wv4l2l1;
    w[4] = Wv5l2l1;
    w[5] = Wv6l2l1;
    w[6] = Wv7l2l1;
    w[7] = Wv8l2l1;
}

/* End of code generation (find_weighting_trilin.cpp) */

/*
               *
               */
void ADCP_measurement_model::eml_sort(const real_T x_data[5], const int32_T x_size[1], real_T
y_data[5], int32_T y_size[1], int32_T idx_data[5], int32_T
idx_size[1])
{
    int32_T dim;
    int32_T vlen;
    real_T vwork_data[5];
    int32_T vstride;
    int32_T k;
    int32_T i1;
    int32_T j;
    int32_T iidx_data[5];
    int32_T i2;
    int32_T idx0_data[5];
    int32_T b_j;
    int32_T pEnd;
    int32_T p;
    int32_T q;
    int32_T qEnd;
    int32_T kEnd;
    dim = 2;
    if (x_size[0] != 1) {
        dim = 1;
    }

    if (dim <= 1) {
        vlen = x_size[0];
    } else {
        vlen = 1;
    }

    y_size[0] = x_size[0];
    idx_size[0] = (int8_T)x_size[0];
    vstride = 1;
    k = 1;
    while (k <= dim - 1) {
        vstride *= x_size[0];
        k = 2;
    }

    i1 = -1;
    for (j = 1; j <= vstride; j++) {
        i1++;
        dim = i1;
        for (k = 0; k < vlen; k++) {
            vwork_data[k] = x_data[dim];
            dim += vstride;
        }

        for (k = 1; k <= (int8_T)vlen; k++) {
            iidx_data[k - 1] = k;
        }

        for (k = 1; k <= (int8_T)vlen - 1; k += 2) {
            if (vwork_data[k - 1] <= vwork_data[k]) {
            } else {
                iidx_data[k - 1] = k + 1;
                iidx_data[k] = k;
            }
        }

        dim = (int8_T)vlen;
        for (i2 = 0; i2 < dim; i2++) {
            idx0_data[i2] = 1;
        }

        dim = 2;
        while (dim < (int8_T)vlen) {
            i2 = dim << 1;
            b_j = 1;
            for (pEnd = 1 + dim; pEnd < (int8_T)vlen + 1; pEnd = qEnd + dim) {
                p = b_j;
                q = pEnd;
                qEnd = b_j + i2;
                if (qEnd > (int8_T)vlen + 1) {
                    qEnd = (int8_T)vlen + 1;
                }

                k = 0;
                kEnd = qEnd - b_j;
                while (k + 1 <= kEnd) {
                    if (vwork_data[iidx_data[p - 1] - 1] <= vwork_data[iidx_data[q - 1] -
                            1]) {
                        idx0_data[k] = iidx_data[p - 1];
                        p++;
                        if (p == pEnd) {
                            while (q < qEnd) {
                                k++;
                                idx0_data[k] = iidx_data[q - 1];
                                q++;
                            }
                        }
                    } else {
                        idx0_data[k] = iidx_data[q - 1];
                        q++;
                        if (q == qEnd) {
                            while (p < pEnd) {
                                k++;
                                idx0_data[k] = iidx_data[p - 1];
                                p++;
                            }
                        }
                    }

                    k++;
                }

                for (k = 0; k + 1 <= kEnd; k++) {
                    iidx_data[(b_j + k) - 1] = idx0_data[k];
                }

                b_j = qEnd;
            }

            dim = i2;
        }

        dim = i1;
        for (k = 0; k < vlen; k++) {
            y_data[dim] = vwork_data[iidx_data[k] - 1];
            idx_data[dim] = iidx_data[k];
            dim += vstride;
        }
    }
}

/*
               * function [current_weightings,grid_loc_row,grid_loc_col,grid_loc_E,failflag] = findInterceptsandWeightings(position,...
               *     euler,depth_cell_size,cell_start,cell_end,blank_d,hori_res,beam_pitch,beam_yaw,P_init,vert_grid_size)
               */

void ADCP_measurement_model::findInterceptsandWeightings(void)
{
    real_T sin_roll;
    real_T cos_roll;
    real_T sin_pitch;
    real_T cos_pitch;
    real_T sin_yaw;
    real_T cos_yaw;
    real_T CBN[9];
    real_T dv0[3];
    int32_T i0;
    real_T ff[3];
    int32_T k;
    real_T overhang;
    int32_T n;
    real_T grid[21];
    real_T w[56];
    int32_T cell;
    boolean_T exitg1;
    real_T y_top;
    real_T y_bottom;
    real_T x_top;
    real_T x_bottom;
    real_T E_top;
    real_T E1;
    real_T i_num;
    real_T intercepts2[21];
    real_T b_x_top[6];
    int32_T intercepts_size_idx_0;
    real_T intercepts_data[21];
    real_T intercepts2_data[5];
    int32_T w_num;
    int32_T intercepts2_size[1];
    int32_T iidx_size[1];
    int32_T iidx_data[5];
    int32_T tmp_size[1];
    real_T tmp_data[5];
    boolean_T exitg2;
    real_T b_grid[3];
    static const int8_T P_init[3] = { 0,0,0 };

    /*  function current_readings = current_readings(current_depth_tops,current_vel,current_bias,depth,range,depth_cell_size) */
    /*  no angular information yet */
    /* 'findInterceptsandWeightings:6' depth = position(3); */
    /* 'findInterceptsandWeightings:7' hori_pos_x = position(1); */
    /* 'findInterceptsandWeightings:8' hori_pos_y = position(2); */
    /* 'findInterceptsandWeightings:9' failflag = 0; */
    failflag = false;

    /* 'findInterceptsandWeightings:11' roll = euler(1); */
    /* 'findInterceptsandWeightings:12' pitch = euler(2); */
    /* 'findInterceptsandWeightings:13' h = euler(3); */
    /* 'findInterceptsandWeightings:15' f = 1; */
    /*  4dof case has no roll or pitch, depth is fixed, beam direction is longer */
    /* 'findInterceptsandWeightings:16' f_N = cos(beam_yaw)*tan(beam_pitch); */
    /*  horizontal correction */
    /* 'findInterceptsandWeightings:17' f_E = sin(beam_yaw)*tan(beam_pitch); */
    /*  horizontal correction */
    /* 'findInterceptsandWeightings:19' ff = b2n_euler([roll pitch h])*[f_N f_E f]'; */
    /*  e2cbn Evaluates the direction cosine matrix by using the Euler angles. */
    /*  The direction Cosine matrix is defined by three successive */
    /*  rotations roll, pitch and yaw which have to occur in that order. */
    /*  This matrix is used to relate vectors in the body frame to the  */
    /*  navigation frame (NED).  */
    /*  */
    /* 	The angles are in radians. */
    /*   */
    /* 'b2n_euler:12' roll = euler(1); */
    /* 'b2n_euler:13' pitch = euler(2); */
    /* 'b2n_euler:14' yaw = euler(3); */
    /* 'b2n_euler:16' sin_roll = sin(roll); */
    sin_roll = sin(euler[0]);

    /* 'b2n_euler:17' cos_roll = cos(roll); */
    cos_roll = cos(euler[0]);

    /* 'b2n_euler:19' sin_pitch = sin(pitch); */
    sin_pitch = sin(euler[1]);

    /* 'b2n_euler:20' cos_pitch = cos(pitch); */
    cos_pitch = cos(euler[1]);

    /* 'b2n_euler:22' sin_yaw = sin(yaw); */
    sin_yaw = sin(euler[2]);

    /* 'b2n_euler:23' cos_yaw = cos(yaw); */
    cos_yaw = cos(euler[2]);

    /* 'b2n_euler:25' CBN = zeros(3,3); */
    /* 'b2n_euler:27' CBN(1,1) =  cos_pitch*cos_yaw; */
    CBN[0] = cos_pitch * cos_yaw;

    /* 'b2n_euler:28' CBN(1,2) = -cos_roll*sin_yaw + sin_roll*sin_pitch*cos_yaw; */
    CBN[3] = -cos_roll * sin_yaw + sin_roll * sin_pitch * cos_yaw;

    /* 'b2n_euler:29' CBN(1,3) =  sin_roll*sin_yaw + cos_roll*sin_pitch*cos_yaw; */
    CBN[6] = sin_roll * sin_yaw + cos_roll * sin_pitch * cos_yaw;

    /* 'b2n_euler:30' CBN(2,1) =  cos_pitch*sin_yaw; */
    CBN[1] = cos_pitch * sin_yaw;

    /* 'b2n_euler:31' CBN(2,2) =  cos_roll*cos_yaw + sin_roll*sin_pitch*sin_yaw; */
    CBN[4] = cos_roll * cos_yaw + sin_roll * sin_pitch * sin_yaw;

    /* 'b2n_euler:32' CBN(2,3) = -sin_roll*cos_yaw + cos_roll*sin_pitch*sin_yaw; */
    CBN[7] = -sin_roll * cos_yaw + cos_roll * sin_pitch * sin_yaw;

    /* 'b2n_euler:33' CBN(3,1) = -sin_pitch; */
    CBN[2] = -sin_pitch;

    /* 'b2n_euler:34' CBN(3,2) =  sin_roll*cos_pitch; */
    CBN[5] = sin_roll * cos_pitch;

    /* 'b2n_euler:35' CBN(3,3) =  cos_roll*cos_pitch; */
    CBN[8] = cos_roll * cos_pitch;
    dv0[0] = cos(beam_yaw) * tan(beam_pitch);
    dv0[1] = sin(beam_yaw) * tan(beam_pitch);
    dv0[2] = 1.0;
    for (i0 = 0; i0 < 3; i0++) {
        ff[i0] = 0.0;
        for (k = 0; k < 3; k++) {
            ff[i0] += CBN[i0 + 3 * k] * dv0[k];
        }
    }

    /* 'findInterceptsandWeightings:20' f = ff(3); */
    /* 'findInterceptsandWeightings:21' f_N = ff(1); */
    /* 'findInterceptsandWeightings:22' f_E = ff(2); */
    /* 'findInterceptsandWeightings:24' overhang = (0.3 + sqrt(0.3^2 - 4*1.7*(-0.075)))/3.4 * depth_cell_size * f; */
    overhang = 0.316057843894554 * depth_cell_size * ff[2];

    /*  15% area overlap total for both sides */
    /* 'findInterceptsandWeightings:26' n = 0; */
    n = -1;

    /* 'findInterceptsandWeightings:27' current_weightings = zeros(8,7,30); */
    memset(&current_weightings[0], 0, 1680U * sizeof(real_T));

    /* 'findInterceptsandWeightings:28' grid_loc_row = zeros(7,30); */
    for (i0 = 0; i0 < 210; i0++) {
        grid_loc_row[i0] = 0.0;

        /* 'findInterceptsandWeightings:29' grid_loc_col = zeros(7,30); */
        grid_loc_col[i0] = 0.0;

        /* 'findInterceptsandWeightings:30' grid_loc_E = zeros(7,30); */
        grid_loc_E[i0] = 0.0;
    }

    /* 'findInterceptsandWeightings:31' grid = zeros(7,3); */
    memset(&grid[0], 0, 21U * sizeof(real_T));

    /* 'findInterceptsandWeightings:32' w = zeros(8,7); */
    memset(&w[0], 0, 56U * sizeof(real_T));

    /* 'findInterceptsandWeightings:34' for cell = cell_start:cell_end */
    cell = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (cell <= (int32_T)(cell_end + (1.0 - cell_start))
                                 - 1)) {
        cos_roll = cell_start + (real_T)cell;

        /* 'findInterceptsandWeightings:36' n = n+1; */
        n++;

        /* 'findInterceptsandWeightings:38' cell_top = depth + f *(blank_d + (cell-1) * depth_cell_size); */
        sin_roll = position[2] + ff[2] * (blank_d + (cos_roll - 1.0) *
                                          depth_cell_size);

        /* 'findInterceptsandWeightings:40' y_top = cell_top + f*depth_cell_size - f*overhang; */
        y_top = (sin_roll + ff[2] * depth_cell_size) - ff[2] * overhang;

        /* 'findInterceptsandWeightings:41' y_bottom = cell_top + 2*f*depth_cell_size + f*overhang; */
        y_bottom = (sin_roll + 2.0 * ff[2] * depth_cell_size) + ff[2] * overhang;

        /* 'findInterceptsandWeightings:43' x_top = hori_pos_x + f_N*(blank_d + (cell-1) * depth_cell_size + depth_cell_size - overhang); */
        x_top = position[0] + ff[0] * (((blank_d + (cos_roll - 1.0) *
                                         depth_cell_size) + depth_cell_size) - overhang);

        /* 'findInterceptsandWeightings:44' x_bottom = hori_pos_x + f_N*(blank_d + (cell-1) * depth_cell_size + 2*depth_cell_size + overhang); */
        x_bottom = position[0] + ff[0] * (((blank_d + (cos_roll - 1.0) *
                                            depth_cell_size) + 2.0 * depth_cell_size) + overhang);

        /* 'findInterceptsandWeightings:46' E_top = hori_pos_y + f_E*(blank_d + (cell-1) * depth_cell_size + depth_cell_size - overhang); */
        E_top = position[1] + ff[1] * (((blank_d + (cos_roll - 1.0) *
                                         depth_cell_size) + depth_cell_size) - overhang);

        /* 'findInterceptsandWeightings:47' E_bottom = hori_pos_y + f_E*(blank_d + (cell-1) * depth_cell_size + 2*depth_cell_size + overhang); */
        sin_roll = position[1] + ff[1] * (((blank_d + (cos_roll - 1.0) *
                                            depth_cell_size) + 2.0 * depth_cell_size) + overhang);

        /* 'findInterceptsandWeightings:50' if x_top < x_bottom */
        if (x_top < x_bottom) {
            /* 'findInterceptsandWeightings:51' x1 = x_top; */
            sin_pitch = x_top;

            /* 'findInterceptsandWeightings:52' x2 = x_bottom; */
            cos_roll = x_bottom;

            /* 'findInterceptsandWeightings:53' y1 = y_top; */
            sin_yaw = y_top;

            /* 'findInterceptsandWeightings:54' y2 = y_bottom; */
            cos_pitch = y_bottom;

            /* 'findInterceptsandWeightings:55' E1 = E_top; */
            E1 = E_top;

            /* 'findInterceptsandWeightings:56' E2 = E_bottom; */
            cos_yaw = sin_roll;
        } else {
            /* 'findInterceptsandWeightings:57' else */
            /* 'findInterceptsandWeightings:58' x2 = x_top; */
            cos_roll = x_top;

            /* 'findInterceptsandWeightings:59' x1 = x_bottom; */
            sin_pitch = x_bottom;

            /* 'findInterceptsandWeightings:60' y2 = y_top; */
            cos_pitch = y_top;

            /* 'findInterceptsandWeightings:61' y1 = y_bottom; */
            sin_yaw = y_bottom;

            /* 'findInterceptsandWeightings:62' E2 = E_top; */
            cos_yaw = E_top;

            /* 'findInterceptsandWeightings:63' E1 = E_bottom; */
            E1 = sin_roll;
        }

        /* 'findInterceptsandWeightings:66' [intercepts2,i_num,failflag] = find_intercepts_3d_4dof_dfki(x1,y1,E1,x2,y2,E2,hori_res,vert_grid_size); */
        find_intercepts_3d_4dof_dfki(sin_pitch, sin_yaw, E1, cos_roll, cos_pitch,
                                     cos_yaw, hori_res, vert_grid_size, intercepts2, &i_num, &failflag);

        /* 'findInterceptsandWeightings:68' if failflag == 1 */
        if (failflag == 1.0) {
            /*          error('x1 > x2') */
            exitg1 = TRUE;
        } else {
            /* 'findInterceptsandWeightings:73' if i_num == 1 */
            if (i_num == 1.0) {
                /* 'findInterceptsandWeightings:74' intercepts = [x_top y_top E_top;x_bottom y_bottom E_bottom]; */
                b_x_top[0] = x_top;
                b_x_top[2] = y_top;
                b_x_top[4] = E_top;
                b_x_top[1] = x_bottom;
                b_x_top[3] = y_bottom;
                b_x_top[5] = sin_roll;
                intercepts_size_idx_0 = 2;
                for (i0 = 0; i0 < 3; i0++) {
                    for (k = 0; k < 2; k++) {
                        intercepts_data[k + (i0 << 1)] = b_x_top[k + (i0 << 1)];
                    }
                }

                /*  easy to localise cell */
            } else {
                /* 'findInterceptsandWeightings:76' else */
                /* 'findInterceptsandWeightings:77' [tmp,indx]=sort(intercepts2(1:i_num-1,2)); */
                w_num = (int32_T)(i_num - 1.0);
                intercepts2_size[0] = (int32_T)(i_num - 1.0);
                for (i0 = 0; i0 < w_num; i0++) {
                    intercepts2_data[i0] = intercepts2[7 + i0];
                }

                eml_sort(intercepts2_data, intercepts2_size, tmp_data, tmp_size,
                         iidx_data, iidx_size);

                /* 'findInterceptsandWeightings:78' intercepts1 = intercepts2(indx,:); */
                /* 'findInterceptsandWeightings:79' intercepts = [x_top y_top E_top;intercepts1;x_bottom y_bottom E_bottom]; */
                intercepts_size_idx_0 = 2 + iidx_size[0];
                intercepts_data[0] = x_top;
                intercepts_data[2 + iidx_size[0]] = y_top;
                intercepts_data[(2 + iidx_size[0]) << 1] = E_top;
                for (i0 = 0; i0 < 3; i0++) {
                    w_num = iidx_size[0];
                    for (k = 0; k < w_num; k++) {
                        intercepts_data[(k + (2 + iidx_size[0]) * i0) + 1] = intercepts2
                                [(iidx_data[k] + 7 * i0) - 1];
                    }
                }

                intercepts_data[1 + iidx_size[0]] = x_bottom;
                intercepts_data[(iidx_size[0] + iidx_size[0]) + 3] = y_bottom;
                intercepts_data[(iidx_size[0] + ((2 + iidx_size[0]) << 1)) + 1] =
                        sin_roll;
            }

            /* 'findInterceptsandWeightings:82' w_num = 1; */
            w_num = 0;

            /* 'findInterceptsandWeightings:83' for k = 1:size(intercepts,1)-1 */
            k = 0;
            exitg2 = FALSE;
            while ((exitg2 == FALSE) && (k <= intercepts_size_idx_0 - 2)) {
                /* 'findInterceptsandWeightings:84' if intercepts(k,2)-(y_top+y_bottom)/2 > intercepts(k+1,2)-(y_top+y_bottom)/2 */
                if (intercepts_data[k + intercepts_size_idx_0] - (y_top + y_bottom) /
                        2.0 > intercepts_data[(k + intercepts_size_idx_0) + 1] - (y_top +
                                                                                  y_bottom) / 2.0) {
                    /*              error('x1 > x2') */
                    /* 'findInterceptsandWeightings:86' failflag = 1; */
                    failflag = 1.0;
                    exitg2 = TRUE;
                } else {
                    /* 'findInterceptsandWeightings:90' x0 = (x_top-floor(intercepts(k,1)/hori_res)*hori_res)/hori_res; */
                    /* 'findInterceptsandWeightings:91' y0 = (E_top-floor(intercepts(k,3)/hori_res)*hori_res)/hori_res; */
                    /* 'findInterceptsandWeightings:92' z0 = (y_top-floor(intercepts(k,2)/vert_grid_size)*vert_grid_size)/vert_grid_size; */
                    /* 'findInterceptsandWeightings:93' xe = (x_bottom-floor(intercepts(k,1)/hori_res)*hori_res)/hori_res; */
                    /* 'findInterceptsandWeightings:94' ye = (E_bottom-floor(intercepts(k,3)/hori_res)*hori_res)/hori_res; */
                    /* 'findInterceptsandWeightings:95' ze = (y_bottom-floor(intercepts(k,2)/vert_grid_size)*vert_grid_size)/vert_grid_size; */
                    /* 'findInterceptsandWeightings:99' l1 = abs((intercepts(k,1)-x_top)/(x_bottom - x_top)); */
                    /* 'findInterceptsandWeightings:100' l2 = abs((intercepts(k+1,1)-x_top)/(x_bottom - x_top)); */
                    /* 'findInterceptsandWeightings:104' w(1:8,w_num) = find_weighting_trilin(xe,ye,ze,x0,y0,z0,l1,l2); */
                    find_weighting_trilin((x_bottom - floor(intercepts_data[k] / hori_res)
                                           * hori_res) / hori_res, (sin_roll - floor(intercepts_data[k +
                                                                    (intercepts_size_idx_0 << 1)] / hori_res) * hori_res) / hori_res,
                            (y_bottom - floor(intercepts_data[k +
                             intercepts_size_idx_0] / vert_grid_size) * vert_grid_size) /
                            vert_grid_size, (x_top - floor(intercepts_data[k]
                                                           / hori_res) * hori_res) / hori_res, (E_top - floor(intercepts_data[k
                                                                                                              + (intercepts_size_idx_0 << 1)] / hori_res) * hori_res) / hori_res,
                            (y_top - floor(intercepts_data[k +
                             intercepts_size_idx_0] / vert_grid_size) * vert_grid_size) /
                            vert_grid_size, fabs((intercepts_data[k] - x_top)
                                                 / (x_bottom - x_top)), fabs((intercepts_data[k + 1] - x_top) /
                            (x_bottom - x_top)), *(real_T (*)[8])&w[w_num << 3]);

                    /* 'findInterceptsandWeightings:106' w_num = w_num + 1; */
                    w_num++;
                    k++;
                }
            }

            /* 'findInterceptsandWeightings:111' for k = 1:size(intercepts,1)-1 */
            for (k = 0; k <= intercepts_size_idx_0 - 2; k++) {
                /* 'findInterceptsandWeightings:112' grid(k,1:3) = find_grid_3d_exp(intercepts(k,1),intercepts(k,2),intercepts(k,3),intercepts(k+1,1),intercepts(k+1,2),intercepts(k+1,3),hori_res,vert_grid_size,hori_res,P_init); */
                /*  finds grid row and column given x1,y1,x2,y2 intercepts or points, and */
                /*  finds the grid it passes through */
                /*  grid = [row column] in grid */
                /* 'find_grid_3d_exp:6' grid = zeros(1,3); */
                /* 'find_grid_3d_exp:8' x_mid = (x1 + x2)/2 - P_init(1); */
                /* 'find_grid_3d_exp:9' y_mid = (y1 + y2)/2; */
                /* 'find_grid_3d_exp:10' z_mid = (z1 + z2)/2  - P_init(2); */
                /* 'find_grid_3d_exp:12' grid(1) = ceil(y_mid/depth_cell_size); */
                b_grid[0] = ceil((intercepts_data[k + intercepts_size_idx_0] +
                                 intercepts_data[(k + intercepts_size_idx_0) + 1]) /
                        2.0 / vert_grid_size);

                /* 'find_grid_3d_exp:13' grid(2) = ceil(x_mid/hori_res); */
                b_grid[1] = ceil(((intercepts_data[k] + intercepts_data[1 + k]) / 2.0 -
                                 P_init[0]) / hori_res);

                /* 'find_grid_3d_exp:14' grid(3) = ceil(z_mid/hori_res2); */
                b_grid[2] = ceil(((intercepts_data[k + (intercepts_size_idx_0 << 1)] +
                                  intercepts_data[(k + (intercepts_size_idx_0 << 1)) +
                                 1]) / 2.0 - P_init[1]) / hori_res);
                for (i0 = 0; i0 < 3; i0++) {
                    grid[k + 7 * i0] = b_grid[i0];
                }
            }

            /* 'findInterceptsandWeightings:115' current_weightings(1:8,1:w_num-1,n) = w(1:8,1:w_num-1); */
            if (1 > w_num) {
                w_num = -1;
            } else {
                w_num--;
            }

            for (i0 = 0; i0 <= w_num; i0++) {
                memcpy(&current_weightings[(i0 << 3) + 56 * n], &w[i0 << 3], sizeof
                        (real_T) << 3);
            }

            /* 'findInterceptsandWeightings:117' grid_loc_row(1:size(intercepts,1)-1,n) = grid(1:size(intercepts,1)-1,1); */
            w_num = intercepts_size_idx_0 - 2;
            for (i0 = 0; i0 <= w_num; i0++) {
                grid_loc_row[i0 + 7 * n] = grid[i0];
                active_states[n]++;
            }

            /* 'findInterceptsandWeightings:118' grid_loc_col(1:size(intercepts,1)-1,n) = grid(1:size(intercepts,1)-1,2); */
            w_num = intercepts_size_idx_0 - 2;
            for (i0 = 0; i0 <= w_num; i0++) {
                grid_loc_col[i0 + 7 * n] = grid[7 + i0];
            }

            /* 'findInterceptsandWeightings:119' grid_loc_E(1:size(intercepts,1)-1,n) = grid(1:size(intercepts,1)-1,3); */
            w_num = intercepts_size_idx_0 - 2;
            for (i0 = 0; i0 <= w_num; i0++) {
                grid_loc_E[i0 + 7 * n] = grid[14 + i0];
            }

            cell++;
        }
    }
}

/* End of code generation (findInterceptsandWeightings.cpp) */

/*
 * function [X_EKF_loc_index,current_measurements,dz_dh,dz_dv,fail_flag] = h_x_ADCP_dfki(cell_index,current_grid_loc,current_weightings_2,grid_loc_row_2,grid_loc_col_2,grid_loc_E_2,...
 *     X_EKF,beam_pitch,beam_yaw,euler,velocity)
 */
double ADCP_measurement_model::h_x_ADCP_dfki(uint16_T cell_index)
{
  real_T sin_roll;
  real_T cos_roll;
  real_T sin_pitch;
  real_T cos_pitch;
  real_T sin_yaw;
  real_T cos_yaw;
  real_T CBN[9];
  boost::array<int,3> array;
  int32_T i;
  real_T ff[3];
  int32_T ix;
  int32_T k;
  boolean_T exitg1;
  real_T vcN;
  real_T vcE;
  real_T vcD;
  int32_T iy;
  boolean_T exitg2;
  static const int8_T grid_matrix_trilin[24] = { 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1 };

  real_T b_array[3];
  real_T sin_beam_yaw;
  real_T sin_heading;
  real_T cos_beam_yaw;
  real_T cos_heading;
  real_T t;
  real_T a;
  real_T b_a;
  real_T J[3];
  real_T c_a;
  real_T d_a;
  real_T e_a;
  real_T b_euler[3];
  real_T b_ff[3];
  real_T b[3];
  unsigned int X_EKF_loc_index;
  double current_measurements; // for the h_x_adcp calc

  /* 'h_x_ADCP_dfki:4' current_weightings=current_weightings_2(1:8,1:7,cell_index); */
  /* 'h_x_ADCP_dfki:5' grid_loc_row=grid_loc_row_2(1:7,cell_index); */
  /* 'h_x_ADCP_dfki:6' grid_loc_col=grid_loc_col_2(1:7,cell_index); */
  /* 'h_x_ADCP_dfki:7' grid_loc_E=grid_loc_E_2(:,cell_index); */
  /* 'h_x_ADCP_dfki:9' current_measurements = zeros(1,1); */
  current_measurements = 0.0;

  /* 'h_x_ADCP_dfki:10' v = velocity; */
  /* 'h_x_ADCP_dfki:11' roll = euler(1); */
  /* 'h_x_ADCP_dfki:12' pitch = euler(2); */
  /* 'h_x_ADCP_dfki:13' h = euler(3); */
  /* 'h_x_ADCP_dfki:15' grid_matrix_trilin = [  0 0 0; */
  /* 'h_x_ADCP_dfki:16'     1 0 0; */
  /* 'h_x_ADCP_dfki:17'     0 1 0; */
  /* 'h_x_ADCP_dfki:18'     0 0 1; */
  /* 'h_x_ADCP_dfki:19'     1 0 1; */
  /* 'h_x_ADCP_dfki:20'     0 1 1; */
  /* 'h_x_ADCP_dfki:21'     1 1 0; */
  /* 'h_x_ADCP_dfki:22'     1 1 1]; */
  /* 'h_x_ADCP_dfki:24' f = 1; */
  /* 'h_x_ADCP_dfki:25' f_N = cos(beam_yaw)*tan(beam_pitch); */
  /*  horizontal correction */
  /* 'h_x_ADCP_dfki:26' f_E = sin(beam_yaw)*tan(beam_pitch); */
  /*  horizontal correction */
  /* 'h_x_ADCP_dfki:28' ff = b2n_euler([roll pitch h])*[f_N f_E f]'; */
  /*  e2cbn Evaluates the direction cosine matrix by using the Euler angles. */
  /*  The direction Cosine matrix is defined by three successive */
  /*  rotations roll, pitch and yaw which have to occur in that order. */
  /*  This matrix is used to relate vectors in the body frame to the  */
  /*  navigation frame (NED).  */
  /*  */
  /* 	The angles are in radians. */
  /*   */
  /* 'b2n_euler:12' roll = euler(1); */
  /* 'b2n_euler:13' pitch = euler(2); */
  /* 'b2n_euler:14' yaw = euler(3); */
  /* 'b2n_euler:16' sin_roll = sin(roll); */
  sin_roll = sin(euler[0]);

  /* 'b2n_euler:17' cos_roll = cos(roll); */
  cos_roll = cos(euler[0]);

  /* 'b2n_euler:19' sin_pitch = sin(pitch); */
  sin_pitch = sin(euler[1]);

  /* 'b2n_euler:20' cos_pitch = cos(pitch); */
  cos_pitch = cos(euler[1]);

  /* 'b2n_euler:22' sin_yaw = sin(yaw); */
  sin_yaw = sin(euler[2]);

  /* 'b2n_euler:23' cos_yaw = cos(yaw); */
  cos_yaw = cos(euler[2]);

  /* 'b2n_euler:25' CBN = zeros(3,3); */
  /* 'b2n_euler:27' CBN(1,1) =  cos_pitch*cos_yaw; */
  CBN[0] = cos_pitch * cos_yaw;

  /* 'b2n_euler:28' CBN(1,2) = -cos_roll*sin_yaw + sin_roll*sin_pitch*cos_yaw; */
  CBN[3] = -cos_roll * sin_yaw + sin_roll * sin_pitch * cos_yaw;

  /* 'b2n_euler:29' CBN(1,3) =  sin_roll*sin_yaw + cos_roll*sin_pitch*cos_yaw; */
  CBN[6] = sin_roll * sin_yaw + cos_roll * sin_pitch * cos_yaw;

  /* 'b2n_euler:30' CBN(2,1) =  cos_pitch*sin_yaw; */
  CBN[1] = cos_pitch * sin_yaw;

  /* 'b2n_euler:31' CBN(2,2) =  cos_roll*cos_yaw + sin_roll*sin_pitch*sin_yaw; */
  CBN[4] = cos_roll * cos_yaw + sin_roll * sin_pitch * sin_yaw;

  /* 'b2n_euler:32' CBN(2,3) = -sin_roll*cos_yaw + cos_roll*sin_pitch*sin_yaw; */
  CBN[7] = -sin_roll * cos_yaw + cos_roll * sin_pitch * sin_yaw;

  /* 'b2n_euler:33' CBN(3,1) = -sin_pitch; */
  CBN[2] = -sin_pitch;

  /* 'b2n_euler:34' CBN(3,2) =  sin_roll*cos_pitch; */
  CBN[5] = sin_roll * cos_pitch;

  /* 'b2n_euler:35' CBN(3,3) =  cos_roll*cos_pitch; */
  CBN[8] = cos_roll * cos_pitch;
  b_array[0] = cos(beam_yaw) * tan(beam_pitch);
  b_array[1] = sin(beam_yaw) * tan(beam_pitch);
  b_array[2] = 1.0;
  for (i = 0; i < 3; i++) {
    ff[i] = 0.0;
    for (ix = 0; ix < 3; ix++) {
      ff[i] += CBN[i + 3 * ix] * b_array[ix];
    }

    /* 'h_x_ADCP_dfki:29' f = ff(3); */
    /* 'h_x_ADCP_dfki:30' f_N = ff(1); */
    /* 'h_x_ADCP_dfki:31' f_E = ff(2); */
    /* 'h_x_ADCP_dfki:32' dz_dh = [0 0 0]; */
    dz_dh[i] = 0.0;
  }

  /* 'h_x_ADCP_dfki:34' fail_flag = 0; */
  failflag = 0.0;

  /* 'h_x_ADCP_dfki:36' k = 1; */
  k = 0;

  /* 'h_x_ADCP_dfki:37' array = zeros(3,1); */
  /* 'h_x_ADCP_dfki:38' X_EKF_loc_index = 0; */
  X_EKF_loc_index = 0.0;

  /* 'h_x_ADCP_dfki:40' while (grid_loc_row(k,1) ~= 0) */
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (grid_loc_row[k + 7 * (cell_index - 1)] != 0.0))
  {
    /* 'h_x_ADCP_dfki:42' vcN = 0; */
    vcN = 0.0;

    /* 'h_x_ADCP_dfki:43' vcE = 0; */
    vcE = 0.0;

    /* 'h_x_ADCP_dfki:44' vcD = 0; */
    vcD = 0.0;

    /* 'h_x_ADCP_dfki:46' for nn = 1:8 */
    iy = 0;
    exitg2 = FALSE;
    while ((exitg2 == FALSE) && (iy < 8)) {
      /* 'h_x_ADCP_dfki:48' array(1) = grid_loc_row(k,1)+grid_matrix_trilin(nn,3); */
      array[0] = grid_loc_row[k + 7 * (cell_index - 1)] + (real_T)
        grid_matrix_trilin[16 + iy];

      /* 'h_x_ADCP_dfki:49' array(2) = grid_loc_col(k,1)+grid_matrix_trilin(nn,1); */
      array[1] = grid_loc_col[k + 7 * (cell_index - 1)] + (real_T)
        grid_matrix_trilin[iy];

      /* 'h_x_ADCP_dfki:50' array(3) = grid_loc_E(k,1)+grid_matrix_trilin(nn,2); */
      array[2] = grid_loc_E[k + 7 * (cell_index - 1)] + (real_T)
        grid_matrix_trilin[8 + iy];

      /* 'h_x_ADCP_dfki:51' X_EKF_loc_index = coder.ceval('Map',array); */

      X_EKF_loc_index = Map[array];

      /*          X_EKF_loc_index = current_grid_loc(grid_loc_row(k,1)+grid_matrix_trilin(nn,3) ... */
      /*              ,grid_loc_col(k,1)+grid_matrix_trilin(nn,1) ... */
      /*              ,grid_loc_E(k,1)+grid_matrix_trilin(nn,2)); % fastest */
      /* 'h_x_ADCP_dfki:58' if X_EKF_loc_index == 0 */
      if (X_EKF_loc_index == 0.0) {
        /*              warning('\nBeam mismatch, rejecting measurements') */
        /* 'h_x_ADCP_dfki:61' fail_flag = 1; */
        failflag = 1.0;
        exitg2 = TRUE;
      } else {
        /* 'h_x_ADCP_dfki:64' else */
        /* 'h_x_ADCP_dfki:66' vcN = vcN + current_weightings(nn,k,1)*(X_EKF(X_EKF_loc_index)-v(1)); */
        i = iy + (k << 3);
        vcN += current_weightings[(i % 8 + (i / 8 << 3)) + 56 * (cell_index -
          1)] * (X_EKF[(int32_T)X_EKF_loc_index - 1] - velocity[0]);

        /* 'h_x_ADCP_dfki:67' vcE = vcE + current_weightings(nn,k,1)*(X_EKF(X_EKF_loc_index+1)-v(2)); */
        i = iy + (k << 3);
        vcE += current_weightings[(i % 8 + (i / 8 << 3)) + 56 * (cell_index -
          1)] * (X_EKF[(int32_T)X_EKF_loc_index] - velocity[1]);

        /* 'h_x_ADCP_dfki:68' vcD = vcD + current_weightings(nn,k,1)*(X_EKF(X_EKF_loc_index+2)-v(3)); */
        i = iy + (k << 3);
        vcD += current_weightings[(i % 8 + (i / 8 << 3)) + 56 * (cell_index -
          1)] * (X_EKF[(int32_T)X_EKF_loc_index + 1] - velocity[2]);

        /* 'h_x_ADCP_dfki:70' vc = [X_EKF(X_EKF_loc_index) X_EKF(X_EKF_loc_index+1) X_EKF(X_EKF_loc_index+2)]'; */
        /* 'h_x_ADCP_dfki:72' dz_dh = dz_dh + current_weightings(nn,k,1)*J_dz_d_att([roll pitch h],beam_yaw,beam_pitch,vc,v); */
        /* 'J_dz_d_att:3' roll = theta_j(1); */
        /* 'J_dz_d_att:4' pitch = theta_j(2); */
        /* 'J_dz_d_att:5' heading = theta_j(3); */
        /* 'J_dz_d_att:7' vx = v(1); */
        /* 'J_dz_d_att:8' vy = v(2); */
        /* 'J_dz_d_att:9' vz = v(3); */
        /* 'J_dz_d_att:11' vxc = vc(1); */
        /* 'J_dz_d_att:12' vyc = vc(2); */
        /* 'J_dz_d_att:13' vzc = vc(3); */
        /* 'J_dz_d_att:15' sin_roll = sin(roll); */
        sin_roll = sin(euler[0]);

        /* 'J_dz_d_att:16' sin_pitch = sin(pitch); */
        sin_pitch = sin(euler[1]);

        /* 'J_dz_d_att:17' sin_beam_pitch = sin(beam_pitch); */
        sin_yaw = sin(beam_pitch);

        /* 'J_dz_d_att:18' sin_beam_yaw = sin(beam_yaw); */
        sin_beam_yaw = sin(beam_yaw);

        /* 'J_dz_d_att:19' sin_heading = sin(heading); */
        sin_heading = sin(euler[2]);

        /* 'J_dz_d_att:21' cos_roll = cos(roll); */
        cos_roll = cos(euler[0]);

        /* 'J_dz_d_att:22' cos_pitch = cos(pitch); */
        cos_pitch = cos(euler[1]);

        /* 'J_dz_d_att:23' cos_beam_pitch = cos(beam_pitch); */
        cos_yaw = cos(beam_pitch);

        /* 'J_dz_d_att:24' cos_beam_yaw = cos(beam_yaw); */
        cos_beam_yaw = cos(beam_yaw);

        /* 'J_dz_d_att:25' cos_heading = cos(heading); */
        cos_heading = cos(euler[2]);

        /* 'J_dz_d_att:27' J = zeros(1,3); */
        /* 'J_dz_d_att:29' J(1,1) = (vy*cos_beam_pitch*cos_heading*cos_roll - vyc*cos_beam_pitch*cos_heading*cos_roll - vx*cos_beam_pitch*cos_roll*sin_heading + vxc*cos_beam_pitch*cos_roll*sin_heading + vz*cos_beam_pitch*cos_pitch*sin_roll - vzc*cos_beam_pitch*cos_pitch*sin_roll - vz*cos_pitch*sin_beam_yaw*sin_beam_pitch*cos_roll + vzc*cos_pitch*sin_beam_yaw*sin_beam_pitch*cos_roll + vx*cos_beam_pitch*cos_heading*sin_pitch*sin_roll - vxc*cos_beam_pitch*cos_heading*sin_pitch*sin_roll + vy*cos_heading*sin_beam_yaw*sin_beam_pitch*sin_roll - vyc*cos_heading*sin_beam_yaw*sin_beam_pitch*sin_roll + vy*cos_beam_pitch*sin_heading*sin_pitch*sin_roll - vyc*cos_beam_pitch*sin_heading*sin_pitch*sin_roll - vx*sin_beam_yaw*sin_beam_pitch*sin_heading*sin_roll + vxc*sin_beam_yaw*sin_beam_pitch*sin_heading*sin_roll - vx*cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll*sin_pitch + vxc*cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll*sin_pitch - vy*sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading*sin_pitch + vyc*sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading*sin_pitch)/(cos_beam_pitch*((cos_beam_pitch*cos_pitch*cos_roll - cos_beam_yaw*sin_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*sin_beam_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_pitch*sin_heading*sin_roll + cos_beam_yaw*cos_heading*cos_pitch*sin_beam_pitch + cos_beam_pitch*cos_heading*cos_roll*sin_pitch - sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading + cos_heading*sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_yaw*cos_pitch*sin_beam_pitch*sin_heading - cos_beam_pitch*cos_heading*sin_roll + cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll + cos_beam_pitch*cos_roll*sin_heading*sin_pitch + sin_beam_yaw*sin_beam_pitch*sin_heading*sin_pitch*sin_roll)^2/cos_beam_pitch^2)^(1/2)); */
        t = (cos_yaw * cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch)
          + cos_pitch * sin_beam_yaw * sin_yaw * sin_roll;
        a = (((cos_yaw * sin_heading * sin_roll + cos_beam_yaw * cos_heading *
               cos_pitch * sin_yaw) + cos_yaw * cos_heading * cos_roll *
              sin_pitch) - sin_beam_yaw * sin_yaw * cos_roll * sin_heading) +
          cos_heading * sin_beam_yaw * sin_yaw * sin_pitch * sin_roll;
        b_a = (((cos_beam_yaw * cos_pitch * sin_yaw * sin_heading - cos_yaw *
                 cos_heading * sin_roll) + cos_heading * sin_beam_yaw * sin_yaw *
                cos_roll) + cos_yaw * cos_roll * sin_heading * sin_pitch) +
          sin_beam_yaw * sin_yaw * sin_heading * sin_pitch * sin_roll;
        J[0] = (((((((((((((((((((velocity[1] * cos_yaw * cos_heading * cos_roll
          - X_EKF[(int32_T)X_EKF_loc_index] * cos_yaw * cos_heading * cos_roll)
          - velocity[0] * cos_yaw * cos_roll * sin_heading) + X_EKF[(int32_T)
          X_EKF_loc_index - 1] * cos_yaw * cos_roll * sin_heading) + velocity[2]
          * cos_yaw * cos_pitch * sin_roll) - X_EKF[(int32_T)X_EKF_loc_index +
                              1] * cos_yaw * cos_pitch * sin_roll) - velocity[2]
                             * cos_pitch * sin_beam_yaw * sin_yaw * cos_roll) +
                            X_EKF[(int32_T)X_EKF_loc_index + 1] * cos_pitch *
                            sin_beam_yaw * sin_yaw * cos_roll) + velocity[0] *
                           cos_yaw * cos_heading * sin_pitch * sin_roll) -
                          X_EKF[(int32_T)X_EKF_loc_index - 1] * cos_yaw *
                          cos_heading * sin_pitch * sin_roll) + velocity[1] *
                         cos_heading * sin_beam_yaw * sin_yaw * sin_roll) -
                        X_EKF[(int32_T)X_EKF_loc_index] * cos_heading *
                        sin_beam_yaw * sin_yaw * sin_roll) + velocity[1] *
                       cos_yaw * sin_heading * sin_pitch * sin_roll) - X_EKF
                      [(int32_T)X_EKF_loc_index] * cos_yaw * sin_heading *
                      sin_pitch * sin_roll) - velocity[0] * sin_beam_yaw *
                     sin_yaw * sin_heading * sin_roll) + X_EKF[(int32_T)
                    X_EKF_loc_index - 1] * sin_beam_yaw * sin_yaw * sin_heading
                    * sin_roll) - velocity[0] * cos_heading * sin_beam_yaw *
                   sin_yaw * cos_roll * sin_pitch) + X_EKF[(int32_T)
                  X_EKF_loc_index - 1] * cos_heading * sin_beam_yaw * sin_yaw *
                  cos_roll * sin_pitch) - velocity[1] * sin_beam_yaw * sin_yaw *
                 cos_roll * sin_heading * sin_pitch) + X_EKF[(int32_T)
                X_EKF_loc_index] * sin_beam_yaw * sin_yaw * cos_roll *
                sin_heading * sin_pitch) / (cos_yaw * sqrt((t * t / (cos_yaw *
          cos_yaw) + a * a / (cos_yaw * cos_yaw)) + b_a * b_a / (cos_yaw *
          cos_yaw)));

        /* 'J_dz_d_att:30' J(1,2) = (abs(cos_beam_pitch)*(vz - vzc)*(cos_beam_yaw*cos_pitch*sin_beam_pitch + cos_beam_pitch*cos_roll*sin_pitch + sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll))/cos_beam_pitch - (cos_heading*(vx - vxc)*(cos_beam_pitch*cos_pitch*cos_roll - cos_beam_yaw*sin_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*sin_beam_pitch*sin_roll))/(cos_beam_pitch*((cos_beam_pitch*cos_pitch*cos_roll - cos_beam_yaw*sin_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*sin_beam_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_pitch*sin_heading*sin_roll + cos_beam_yaw*cos_heading*cos_pitch*sin_beam_pitch + cos_beam_pitch*cos_heading*cos_roll*sin_pitch - sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading + cos_heading*sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_yaw*cos_pitch*sin_beam_pitch*sin_heading - cos_beam_pitch*cos_heading*sin_roll + cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll + cos_beam_pitch*cos_roll*sin_heading*sin_pitch + sin_beam_yaw*sin_beam_pitch*sin_heading*sin_pitch*sin_roll)^2/cos_beam_pitch^2)^(1/2)) - (sin_heading*(vy - vyc)*(cos_beam_pitch*cos_pitch*cos_roll - cos_beam_yaw*sin_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*sin_beam_pitch*sin_roll))/(cos_beam_pitch*((cos_beam_pitch*cos_pitch*cos_roll - cos_beam_yaw*sin_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*sin_beam_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_pitch*sin_heading*sin_roll + cos_beam_yaw*cos_heading*cos_pitch*sin_beam_pitch + cos_beam_pitch*cos_heading*cos_roll*sin_pitch - sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading + cos_heading*sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_yaw*cos_pitch*sin_beam_pitch*sin_heading - cos_beam_pitch*cos_heading*sin_roll + cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll + cos_beam_pitch*cos_roll*sin_heading*sin_pitch + sin_beam_yaw*sin_beam_pitch*sin_heading*sin_pitch*sin_roll)^2/cos_beam_pitch^2)^(1/2)); */
        t = (cos_yaw * cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch)
          + cos_pitch * sin_beam_yaw * sin_yaw * sin_roll;
        a = (((cos_yaw * sin_heading * sin_roll + cos_beam_yaw * cos_heading *
               cos_pitch * sin_yaw) + cos_yaw * cos_heading * cos_roll *
              sin_pitch) - sin_beam_yaw * sin_yaw * cos_roll * sin_heading) +
          cos_heading * sin_beam_yaw * sin_yaw * sin_pitch * sin_roll;
        b_a = (((cos_beam_yaw * cos_pitch * sin_yaw * sin_heading - cos_yaw *
                 cos_heading * sin_roll) + cos_heading * sin_beam_yaw * sin_yaw *
                cos_roll) + cos_yaw * cos_roll * sin_heading * sin_pitch) +
          sin_beam_yaw * sin_yaw * sin_heading * sin_pitch * sin_roll;
        c_a = (cos_yaw * cos_pitch * cos_roll - cos_beam_yaw * sin_yaw *
               sin_pitch) + cos_pitch * sin_beam_yaw * sin_yaw * sin_roll;
        d_a = (((cos_yaw * sin_heading * sin_roll + cos_beam_yaw * cos_heading *
                 cos_pitch * sin_yaw) + cos_yaw * cos_heading * cos_roll *
                sin_pitch) - sin_beam_yaw * sin_yaw * cos_roll * sin_heading) +
          cos_heading * sin_beam_yaw * sin_yaw * sin_pitch * sin_roll;
        e_a = (((cos_beam_yaw * cos_pitch * sin_yaw * sin_heading - cos_yaw *
                 cos_heading * sin_roll) + cos_heading * sin_beam_yaw * sin_yaw *
                cos_roll) + cos_yaw * cos_roll * sin_heading * sin_pitch) +
          sin_beam_yaw * sin_yaw * sin_heading * sin_pitch * sin_roll;
        J[1] = (fabs(cos_yaw) * (velocity[2] - X_EKF[(int32_T)X_EKF_loc_index +
                 1]) * ((cos_beam_yaw * cos_pitch * sin_yaw + cos_yaw * cos_roll
                         * sin_pitch) + sin_beam_yaw * sin_yaw * sin_pitch *
                        sin_roll) / cos_yaw - cos_heading * (velocity[0] -
                 X_EKF[(int32_T)X_EKF_loc_index - 1]) * ((cos_yaw * cos_pitch *
                  cos_roll - cos_beam_yaw * sin_yaw * sin_pitch) + cos_pitch *
                 sin_beam_yaw * sin_yaw * sin_roll) / (cos_yaw * sqrt((t * t /
                   (cos_yaw * cos_yaw) + a * a / (cos_yaw * cos_yaw)) + b_a *
                  b_a / (cos_yaw * cos_yaw)))) - sin_heading * (velocity[1] -
          X_EKF[(int32_T)X_EKF_loc_index]) * ((cos_yaw * cos_pitch * cos_roll -
          cos_beam_yaw * sin_yaw * sin_pitch) + cos_pitch * sin_beam_yaw *
          sin_yaw * sin_roll) / (cos_yaw * sqrt((c_a * c_a / (cos_yaw * cos_yaw)
          + d_a * d_a / (cos_yaw * cos_yaw)) + e_a * e_a / (cos_yaw * cos_yaw)));

        /* 'J_dz_d_att:31' J(1,3) = -(vx*cos_beam_pitch*cos_heading*sin_roll - vxc*cos_beam_pitch*cos_heading*sin_roll + vy*cos_beam_pitch*sin_heading*sin_roll - vyc*cos_beam_pitch*sin_heading*sin_roll + vy*cos_beam_yaw*cos_heading*cos_pitch*sin_beam_pitch - vyc*cos_beam_yaw*cos_heading*cos_pitch*sin_beam_pitch - vx*cos_beam_yaw*cos_pitch*sin_beam_pitch*sin_heading + vxc*cos_beam_yaw*cos_pitch*sin_beam_pitch*sin_heading + vy*cos_beam_pitch*cos_heading*cos_roll*sin_pitch - vyc*cos_beam_pitch*cos_heading*cos_roll*sin_pitch - vx*cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll + vxc*cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll - vx*cos_beam_pitch*cos_roll*sin_heading*sin_pitch + vxc*cos_beam_pitch*cos_roll*sin_heading*sin_pitch - vy*sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading + vyc*sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading + vy*cos_heading*sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll - vyc*cos_heading*sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll - vx*sin_beam_yaw*sin_beam_pitch*sin_heading*sin_pitch*sin_roll + vxc*sin_beam_yaw*sin_beam_pitch*sin_heading*sin_pitch*sin_roll)/(cos_beam_pitch*((cos_beam_pitch*cos_pitch*cos_roll - cos_beam_yaw*sin_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*sin_beam_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_pitch*sin_heading*sin_roll + cos_beam_yaw*cos_heading*cos_pitch*sin_beam_pitch + cos_beam_pitch*cos_heading*cos_roll*sin_pitch - sin_beam_yaw*sin_beam_pitch*cos_roll*sin_heading + cos_heading*sin_beam_yaw*sin_beam_pitch*sin_pitch*sin_roll)^2/cos_beam_pitch^2 + (cos_beam_yaw*cos_pitch*sin_beam_pitch*sin_heading - cos_beam_pitch*cos_heading*sin_roll + cos_heading*sin_beam_yaw*sin_beam_pitch*cos_roll + cos_beam_pitch*cos_roll*sin_heading*sin_pitch + sin_beam_yaw*sin_beam_pitch*sin_heading*sin_pitch*sin_roll)^2/cos_beam_pitch^2)^(1/2)); */
        t = (cos_yaw * cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch)
          + cos_pitch * sin_beam_yaw * sin_yaw * sin_roll;
        a = (((cos_yaw * sin_heading * sin_roll + cos_beam_yaw * cos_heading *
               cos_pitch * sin_yaw) + cos_yaw * cos_heading * cos_roll *
              sin_pitch) - sin_beam_yaw * sin_yaw * cos_roll * sin_heading) +
          cos_heading * sin_beam_yaw * sin_yaw * sin_pitch * sin_roll;
        b_a = (((cos_beam_yaw * cos_pitch * sin_yaw * sin_heading - cos_yaw *
                 cos_heading * sin_roll) + cos_heading * sin_beam_yaw * sin_yaw *
                cos_roll) + cos_yaw * cos_roll * sin_heading * sin_pitch) +
          sin_beam_yaw * sin_yaw * sin_heading * sin_pitch * sin_roll;
        J[2] = -(((((((((((((((((((velocity[0] * cos_yaw * cos_heading *
          sin_roll - X_EKF[(int32_T)X_EKF_loc_index - 1] * cos_yaw *
          cos_heading * sin_roll) + velocity[1] * cos_yaw * sin_heading *
          sin_roll) - X_EKF[(int32_T)X_EKF_loc_index] * cos_yaw * sin_heading *
          sin_roll) + velocity[1] * cos_beam_yaw * cos_heading * cos_pitch *
          sin_yaw) - X_EKF[(int32_T)X_EKF_loc_index] * cos_beam_yaw *
          cos_heading * cos_pitch * sin_yaw) - velocity[0] * cos_beam_yaw *
                              cos_pitch * sin_yaw * sin_heading) + X_EKF
                             [(int32_T)X_EKF_loc_index - 1] * cos_beam_yaw *
                             cos_pitch * sin_yaw * sin_heading) + velocity[1] *
                            cos_yaw * cos_heading * cos_roll * sin_pitch) -
                           X_EKF[(int32_T)X_EKF_loc_index] * cos_yaw *
                           cos_heading * cos_roll * sin_pitch) - velocity[0] *
                          cos_heading * sin_beam_yaw * sin_yaw * cos_roll) +
                         X_EKF[(int32_T)X_EKF_loc_index - 1] * cos_heading *
                         sin_beam_yaw * sin_yaw * cos_roll) - velocity[0] *
                        cos_yaw * cos_roll * sin_heading * sin_pitch) + X_EKF
                       [(int32_T)X_EKF_loc_index - 1] * cos_yaw * cos_roll *
                       sin_heading * sin_pitch) - velocity[1] * sin_beam_yaw *
                      sin_yaw * cos_roll * sin_heading) + X_EKF[(int32_T)
                     X_EKF_loc_index] * sin_beam_yaw * sin_yaw * cos_roll *
                     sin_heading) + velocity[1] * cos_heading * sin_beam_yaw *
                    sin_yaw * sin_pitch * sin_roll) - X_EKF[(int32_T)
                   X_EKF_loc_index] * cos_heading * sin_beam_yaw * sin_yaw *
                   sin_pitch * sin_roll) - velocity[0] * sin_beam_yaw * sin_yaw *
                  sin_heading * sin_pitch * sin_roll) + X_EKF[(int32_T)
                 X_EKF_loc_index - 1] * sin_beam_yaw * sin_yaw * sin_heading *
                 sin_pitch * sin_roll) / (cos_yaw * sqrt((t * t / (cos_yaw *
          cos_yaw) + a * a / (cos_yaw * cos_yaw)) + b_a * b_a / (cos_yaw *
          cos_yaw)));

        /*  offset = 1e-9; */
        /*   */
        /*  y = feval(@dz_ADCP,theta_j,beam_yaw,beam_pitch,vc,v); */
        /*   */
        /*  lenx = length(theta_j); */
        /*  leny = length(y); */
        /*   */
        /*  J = zeros(leny,lenx); */
        /*   */
        /*  xt2 = theta_j; */
        /*   */
        /*  for i=1:lenx */
        /*      xt2(i) = theta_j(i) + offset; */
        /*      yt = feval(@dz_ADCP,xt2,beam_yaw,beam_pitch,vc,v); */
        /*      dy = feval(@default_dmodel, yt, y); */
        /*      J(:,i) = dy/offset; */
        /*      xt2(i) = theta_j(i); */
        /*  end */
        i = iy + (k << 3);
        for (ix = 0; ix < 3; ix++) {
          dz_dh[ix] += current_weightings[(i % 8 + (i / 8 << 3)) + 56 *
            (cell_index - 1)] * J[ix];
        }

        iy++;
      }
    }

    /* 'h_x_ADCP_dfki:78' current_measurements(1,1) = current_measurements(1,1) + dot([f_N f_E f]/norm([f_N f_E f]),[vcN vcE vcD]); */
    b_euler[0] = ff[0];
    b_euler[1] = ff[1];
    b_euler[2] = ff[2];
    a = 0.0;
    sin_yaw = 2.2250738585072014E-308;
    for (i = 0; i < 3; i++) {
      cos_yaw = fabs(b_euler[i]);
      if (cos_yaw > sin_yaw) {
        t = sin_yaw / cos_yaw;
        a = 1.0 + a * t * t;
        sin_yaw = cos_yaw;
      } else {
        t = cos_yaw / sin_yaw;
        a += t * t;
      }
    }

    a = sin_yaw * sqrt(a);
    b_ff[0] = ff[0];
    b_ff[1] = ff[1];
    b_ff[2] = ff[2];
    for (i = 0; i < 3; i++) {
      b_euler[i] = b_ff[i] / a;
    }

    b[0] = vcN;
    b[1] = vcE;
    b[2] = vcD;
    sin_yaw = 0.0;
    ix = 0;
    iy = 0;
    for (i = 0; i < 3; i++) {
      sin_yaw += b_euler[ix] * b[iy];
      ix++;
      iy++;
    }

    current_measurements += sin_yaw;

    /* 'h_x_ADCP_dfki:79' k = k + 1; */
    k++;

    /* 'h_x_ADCP_dfki:81' if k > size(grid_loc_row,1) */
    if (k + 1 > 7) {
      exitg1 = TRUE;
    }
  }

  /* 'h_x_ADCP_dfki:87' dz_dv = J_dz_dv([roll pitch h],beam_yaw,beam_pitch); */
  /* 'J_dz_dv:3' roll = theta_j(1); */
  /* 'J_dz_dv:4' pitch = theta_j(2); */
  /* 'J_dz_dv:5' heading = theta_j(3); */
  /* 'J_dz_dv:7' sin_roll = sin(roll); */
  sin_roll = sin(euler[0]);

  /* 'J_dz_dv:8' sin_pitch = sin(pitch); */
  sin_pitch = sin(euler[1]);

  /* 'J_dz_dv:9' sin_beam_pitch = sin(beam_pitch); */
  /* 'J_dz_dv:10' sin_beam_yaw = sin(beam_yaw); */
  sin_beam_yaw = sin(beam_yaw);

  /* 'J_dz_dv:11' sin_heading = sin(heading); */
  sin_heading = sin(euler[2]);

  /* 'J_dz_dv:13' cos_roll = cos(roll); */
  cos_roll = cos(euler[0]);

  /* 'J_dz_dv:14' cos_pitch = cos(pitch); */
  cos_pitch = cos(euler[1]);

  /* 'J_dz_dv:15' cos_beam_pitch = cos(beam_pitch); */
  /* 'J_dz_dv:16' cos_beam_yaw = cos(beam_yaw); */
  cos_beam_yaw = cos(beam_yaw);

  /* 'J_dz_dv:17' cos_heading = cos(heading); */
  cos_heading = cos(euler[2]);

  /* 'J_dz_dv:19' tan_beam_pitch = sin_beam_pitch/cos_beam_pitch; */
  sin_yaw = sin(beam_pitch) / cos(beam_pitch);

  /* 'J_dz_dv:21' J=zeros(1,3); */
  /* 'J_dz_dv:23' J(1,1) = -(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)/(abs(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)^2 + abs(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)^2 + abs(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)^2)^(1/2); */
  a = fabs((cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch) +
           cos_pitch * sin_beam_yaw * sin_yaw * sin_roll);
  cos_yaw = fabs(((sin_heading * sin_roll - sin_beam_yaw * sin_yaw * (cos_roll *
    sin_heading - cos_heading * sin_pitch * sin_roll)) + cos_heading * cos_roll *
                  sin_pitch) + cos_beam_yaw * cos_heading * cos_pitch * sin_yaw);
  t = fabs(((sin_beam_yaw * sin_yaw * (cos_heading * cos_roll + sin_heading *
              sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll *
            sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch * sin_yaw *
           sin_heading);
  dz_dv[0] = -(((sin_heading * sin_roll - sin_beam_yaw * sin_yaw * (cos_roll *
    sin_heading - cos_heading * sin_pitch * sin_roll)) + cos_heading * cos_roll *
                sin_pitch) + cos_beam_yaw * cos_heading * cos_pitch * sin_yaw) /
    sqrt((a * a + cos_yaw * cos_yaw) + t * t);

  /* 'J_dz_dv:24' J(1,2) = -(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)/(abs(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)^2 + abs(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)^2 + abs(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)^2)^(1/2); */
  a = fabs((cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch) +
           cos_pitch * sin_beam_yaw * sin_yaw * sin_roll);
  cos_yaw = fabs(((sin_heading * sin_roll - sin_beam_yaw * sin_yaw * (cos_roll *
    sin_heading - cos_heading * sin_pitch * sin_roll)) + cos_heading * cos_roll *
                  sin_pitch) + cos_beam_yaw * cos_heading * cos_pitch * sin_yaw);
  t = fabs(((sin_beam_yaw * sin_yaw * (cos_heading * cos_roll + sin_heading *
              sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll *
            sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch * sin_yaw *
           sin_heading);
  dz_dv[1] = -(((sin_beam_yaw * sin_yaw * (cos_heading * cos_roll + sin_heading *
    sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll * sin_heading *
                sin_pitch) + cos_beam_yaw * cos_pitch * sin_yaw * sin_heading) /
    sqrt((a * a + cos_yaw * cos_yaw) + t * t);

  /* 'J_dz_dv:25' J(1,3) = -(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)/(abs(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)^2 + abs(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)^2 + abs(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)^2)^(1/2); */
  a = fabs((cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch) +
           cos_pitch * sin_beam_yaw * sin_yaw * sin_roll);
  cos_yaw = fabs(((sin_heading * sin_roll - sin_beam_yaw * sin_yaw * (cos_roll *
    sin_heading - cos_heading * sin_pitch * sin_roll)) + cos_heading * cos_roll *
                  sin_pitch) + cos_beam_yaw * cos_heading * cos_pitch * sin_yaw);
  t = fabs(((sin_beam_yaw * sin_yaw * (cos_heading * cos_roll + sin_heading *
              sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll *
            sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch * sin_yaw *
           sin_heading);
  dz_dv[2] = -((cos_pitch * cos_roll - cos_beam_yaw * sin_yaw * sin_pitch) +
               cos_pitch * sin_beam_yaw * sin_yaw * sin_roll) / sqrt((a * a +
    cos_yaw * cos_yaw) + t * t);

  /*  offset = 1e-9; */
  /*   */
  /*  y = feval(@dz_ADCP,theta_j,beam_yaw,beam_pitch,vc,v); */
  /*   */
  /*  lenx = length(v); */
  /*  leny = length(y); */
  /*   */
  /*  J = zeros(leny,lenx); */
  /*   */
  /*  xt2 = v; */
  /*   */
  /*  for i=1:lenx */
  /*      xt2(i) = v(i) + offset; */
  /*      yt = feval(@dz_ADCP,theta_j,beam_yaw,beam_pitch,vc,xt2); */
  /*      dy = feval(@default_dmodel, yt, y); */
  /*      J(:,i) = dy/offset; */
  /*      xt2(i) = v(i); */
  /*  end */
  return current_measurements;
}

/* End of code generation (h_x_ADCP_dfki.cpp) */

/*
 * function [cell_loc_index,J_vc] = calculateWaterCurrentVelocityJacobian(cell_index,current_grid_loc,current_weightings,grid_loc_row,grid_loc_col,grid_loc_E,euler,beam_pitch,beam_yaw)
 */
void ADCP_measurement_model::calculateWaterCurrentVelocityJacobian(uint16_T cell_index)
{
  int32_T kk;
  boolean_T exitg1;
  int32_T nn_cell_t;
  static const int8_T grid_matrix_trilin[24] = { 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1 };

  boost::array<int,3> array;
  real_T b_array[3];
  int32_T i;
  real_T sin_roll;
  real_T sin_pitch;
  real_T sin_beam_yaw;
  real_T sin_heading;
  real_T cos_roll;
  real_T cos_pitch;
  real_T cos_beam_yaw;
  real_T cos_heading;
  real_T tan_beam_pitch;
  real_T y;
  real_T b_y;
  real_T c_y;
  real_T J[3];
  unsigned int cell_loc_index = 0;

  /* 'calculateWaterCurrentVelocityJacobian:3' roll = euler(1); */
  /* 'calculateWaterCurrentVelocityJacobian:4' pitch = euler (2); */
  /* 'calculateWaterCurrentVelocityJacobian:5' yaw = euler(3); */
  /* 'calculateWaterCurrentVelocityJacobian:7' grid_matrix_trilin = [  0 0 0; */
  /* 'calculateWaterCurrentVelocityJacobian:8'     1 0 0; */
  /* 'calculateWaterCurrentVelocityJacobian:9'     0 1 0; */
  /* 'calculateWaterCurrentVelocityJacobian:10'     0 0 1; */
  /* 'calculateWaterCurrentVelocityJacobian:11'     1 0 1; */
  /* 'calculateWaterCurrentVelocityJacobian:12'     0 1 1; */
  /* 'calculateWaterCurrentVelocityJacobian:13'     1 1 0; */
  /* 'calculateWaterCurrentVelocityJacobian:14'     1 1 1]; */
  /* 'calculateWaterCurrentVelocityJacobian:16' kk = 1; */
  kk = 0;

  /* 'calculateWaterCurrentVelocityJacobian:18' J_vc = zeros(1,200); */
  memset(&J_vc[0], 0, 200U * sizeof(real_T));

  /* 'calculateWaterCurrentVelocityJacobian:20' cell_loc_index = 0; */
  cell_loc_index = 0;

  /* 'calculateWaterCurrentVelocityJacobian:21' array = zeros(3,1); */
  /* 'calculateWaterCurrentVelocityJacobian:25' while grid_loc_row(kk,cell_index) ~= 0 */
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (grid_loc_row[kk + 7 * (cell_index - 1)] != 0.0))
  {
    /*  for each of the grids involved */
    /* 'calculateWaterCurrentVelocityJacobian:27' for nn_cell_t = 1:8 */
    for (nn_cell_t = 0; nn_cell_t < 8; nn_cell_t++) {
      /* 'calculateWaterCurrentVelocityJacobian:29' array(1) = grid_loc_row(kk,cell_index)+grid_matrix_trilin(nn_cell_t,3); */
      array[0] = grid_loc_row[kk + 7 * (cell_index - 1)] + (real_T)
        grid_matrix_trilin[16 + nn_cell_t];

      /* 'calculateWaterCurrentVelocityJacobian:30' array(2) = grid_loc_col(kk,cell_index)+grid_matrix_trilin(nn_cell_t,1); */
      array[1] = grid_loc_col[kk + 7 * (cell_index - 1)] + (real_T)
        grid_matrix_trilin[nn_cell_t];

      /* 'calculateWaterCurrentVelocityJacobian:31' array(3) = grid_loc_E(kk,cell_index)+grid_matrix_trilin(nn_cell_t,2); */
      array[2] = grid_loc_E[kk + 7 * (cell_index - 1)] + (real_T)
        grid_matrix_trilin[8 + nn_cell_t];

      cell_loc_index = Map[array];

      /*          cell_loc_index = current_grid_loc(grid_loc_row(kk,cell_index)+grid_matrix_trilin(nn_cell_t,3)... */
      /*              ,grid_loc_col(kk,cell_index)+grid_matrix_trilin(nn_cell_t,1)... */
      /*              ,grid_loc_E(kk,cell_index)+grid_matrix_trilin(nn_cell_t,2)); */
      /* 'calculateWaterCurrentVelocityJacobian:37' J_vc(1,cell_loc_index:cell_loc_index+2) = J_vc(1,cell_loc_index:cell_loc_index+2) + J_Zadcp_vc(roll,pitch,yaw,beam_yaw,beam_pitch,current_weightings(nn_cell_t,kk,cell_index)); */
      /* 'J_Zadcp_vc:3' sin_roll = sin(roll); */
      sin_roll = sin(euler[0]);

      /* 'J_Zadcp_vc:4' sin_pitch = sin(pitch); */
      sin_pitch = sin(euler[1]);

      /* 'J_Zadcp_vc:5' sin_beam_pitch = sin(beam_pitch); */
      /* 'J_Zadcp_vc:6' sin_beam_yaw = sin(beam_yaw); */
      sin_beam_yaw = sin(beam_yaw);

      /* 'J_Zadcp_vc:7' sin_heading = sin(heading); */
      sin_heading = sin(euler[2]);

      /* 'J_Zadcp_vc:9' cos_roll = cos(roll); */
      cos_roll = cos(euler[0]);

      /* 'J_Zadcp_vc:10' cos_pitch = cos(pitch); */
      cos_pitch = cos(euler[1]);

      /* 'J_Zadcp_vc:11' cos_beam_pitch = cos(beam_pitch); */
      /* 'J_Zadcp_vc:12' cos_beam_yaw = cos(beam_yaw); */
      cos_beam_yaw = cos(beam_yaw);

      /* 'J_Zadcp_vc:13' cos_heading = cos(heading); */
      cos_heading = cos(euler[2]);

      /*  tan_roll = sin_roll/cos_roll; */
      /*  tan_pitch = sin_pitch/cos_pitch; */
      /* 'J_Zadcp_vc:17' tan_beam_pitch = sin_beam_pitch/cos_beam_pitch; */
      tan_beam_pitch = sin(beam_pitch) / cos(beam_pitch);

      /*  tan_beam_yaw = sin_beam_yaw/cos_beam_yaw; */
      /*  tan_heading = sin_heading/cos_heading; */
      /* 'J_Zadcp_vc:21' J = zeros(1,3); */
      /* 'J_Zadcp_vc:23' J(1,1) = c1*((sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)/(abs(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)^2 + abs(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)^2 + abs(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)^2)^(1/2)); */
      y = fabs((cos_pitch * cos_roll - cos_beam_yaw * tan_beam_pitch * sin_pitch)
               + cos_pitch * sin_beam_yaw * tan_beam_pitch * sin_roll);
      b_y = fabs(((sin_heading * sin_roll - sin_beam_yaw * tan_beam_pitch *
                   (cos_roll * sin_heading - cos_heading * sin_pitch * sin_roll))
                  + cos_heading * cos_roll * sin_pitch) + cos_beam_yaw *
                 cos_heading * cos_pitch * tan_beam_pitch);
      c_y = fabs(((sin_beam_yaw * tan_beam_pitch * (cos_heading * cos_roll +
        sin_heading * sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll
                  * sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch *
                 tan_beam_pitch * sin_heading);
      J[0] = current_weightings[(nn_cell_t + (kk << 3)) + 56 * (cell_index - 1)]
        * ((((sin_heading * sin_roll - sin_beam_yaw * tan_beam_pitch * (cos_roll
               * sin_heading - cos_heading * sin_pitch * sin_roll)) +
             cos_heading * cos_roll * sin_pitch) + cos_beam_yaw * cos_heading *
            cos_pitch * tan_beam_pitch) / sqrt((y * y + b_y * b_y) + c_y * c_y));

      /* 'J_Zadcp_vc:24' J(1,2) = c1*((sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)/(abs(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)^2 + abs(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)^2 + abs(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)^2)^(1/2)); */
      y = fabs((cos_pitch * cos_roll - cos_beam_yaw * tan_beam_pitch * sin_pitch)
               + cos_pitch * sin_beam_yaw * tan_beam_pitch * sin_roll);
      b_y = fabs(((sin_heading * sin_roll - sin_beam_yaw * tan_beam_pitch *
                   (cos_roll * sin_heading - cos_heading * sin_pitch * sin_roll))
                  + cos_heading * cos_roll * sin_pitch) + cos_beam_yaw *
                 cos_heading * cos_pitch * tan_beam_pitch);
      c_y = fabs(((sin_beam_yaw * tan_beam_pitch * (cos_heading * cos_roll +
        sin_heading * sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll
                  * sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch *
                 tan_beam_pitch * sin_heading);
      J[1] = current_weightings[(nn_cell_t + (kk << 3)) + 56 * (cell_index - 1)]
        * ((((sin_beam_yaw * tan_beam_pitch * (cos_heading * cos_roll +
               sin_heading * sin_pitch * sin_roll) - cos_heading * sin_roll) +
             cos_roll * sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch *
            tan_beam_pitch * sin_heading) / sqrt((y * y + b_y * b_y) + c_y * c_y));

      /* 'J_Zadcp_vc:25' J(1,3) = c1*((cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)/(abs(cos_pitch*cos_roll - cos_beam_yaw*tan_beam_pitch*sin_pitch + cos_pitch*sin_beam_yaw*tan_beam_pitch*sin_roll)^2 + abs(sin_heading*sin_roll - sin_beam_yaw*tan_beam_pitch*(cos_roll*sin_heading - cos_heading*sin_pitch*sin_roll) + cos_heading*cos_roll*sin_pitch + cos_beam_yaw*cos_heading*cos_pitch*tan_beam_pitch)^2 + abs(sin_beam_yaw*tan_beam_pitch*(cos_heading*cos_roll + sin_heading*sin_pitch*sin_roll) - cos_heading*sin_roll + cos_roll*sin_heading*sin_pitch + cos_beam_yaw*cos_pitch*tan_beam_pitch*sin_heading)^2)^(1/2)); */
      y = fabs((cos_pitch * cos_roll - cos_beam_yaw * tan_beam_pitch * sin_pitch)
               + cos_pitch * sin_beam_yaw * tan_beam_pitch * sin_roll);
      b_y = fabs(((sin_heading * sin_roll - sin_beam_yaw * tan_beam_pitch *
                   (cos_roll * sin_heading - cos_heading * sin_pitch * sin_roll))
                  + cos_heading * cos_roll * sin_pitch) + cos_beam_yaw *
                 cos_heading * cos_pitch * tan_beam_pitch);
      c_y = fabs(((sin_beam_yaw * tan_beam_pitch * (cos_heading * cos_roll +
        sin_heading * sin_pitch * sin_roll) - cos_heading * sin_roll) + cos_roll
                  * sin_heading * sin_pitch) + cos_beam_yaw * cos_pitch *
                 tan_beam_pitch * sin_heading);
      J[2] = current_weightings[(nn_cell_t + (kk << 3)) + 56 * (cell_index - 1)]
        * (((cos_pitch * cos_roll - cos_beam_yaw * tan_beam_pitch * sin_pitch) +
            cos_pitch * sin_beam_yaw * tan_beam_pitch * sin_roll) / sqrt((y * y
             + b_y * b_y) + c_y * c_y));
      for (i = 0; i < 3; i++) {
        b_array[i] = J_vc[(int32_T)(cell_loc_index + (real_T)i) - 1] + J[i];
      }

      for (i = 0; i < 3; i++) {
        J_vc[(int32_T)(cell_loc_index + (real_T)i) - 1] = b_array[i];
      }
    }

    /* 'calculateWaterCurrentVelocityJacobian:42' kk = kk + 1; */
    kk++;

    /* 'calculateWaterCurrentVelocityJacobian:44' if kk > size(grid_loc_row,1) */
    if (kk + 1 > 7) {
      exitg1 = TRUE;
    }
  }
}

/* End of code generation (calculateWaterCurrentVelocityJacobian.cpp) */

void ADCP_measurement_model::test(void)
{
    std::cout << "Hello World!" << std::endl;
    //    std::cout << (&X_EKF)[10] << std::endl;
    std::cout << X_EKF[199] << std::endl;
    std::cout << "Position [2]: ";
    std::cout << position[2] << std::endl;

}

