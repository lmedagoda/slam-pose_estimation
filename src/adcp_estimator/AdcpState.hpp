#ifndef _ADCP_STATE_HPP_
#define _ADCP_STATE_HPP_

/** MTK library **/
#include <mtk/src/SubManifold.hpp>
#include <mtk/types/SOn.hpp>
#include <mtk/types/vect.hpp>
#include <mtk/build_manifold.hpp>
#include <mtk/startIdx.hpp>
#include <boost/concept_check.hpp>
#include <pose_estimation/ManifoldHelper.hpp>
#include <Eigen/Core>

namespace pose_estimation
{

typedef MTK::SO3<double> RotationType;
typedef RotationType::vect_type TranslationType;
typedef RotationType::vect_type VelocityType;
typedef RotationType::vect_type BiasType;
typedef MTK::vect<4, double> ADCPBiasType;
typedef MTK::vect<60, double> VCurrentType;

MTK_BUILD_MANIFOLD(AdcpState,
   ((RotationType, orientation))
   ((VelocityType, velocity)) // navigation/target frame velocity
   ((TranslationType, position)) // navigation/target frame velocity
   ((BiasType, bias_gyro))
   ((BiasType, bias_acc))
   ((ADCPBiasType, bias_adcp))
   ((VCurrentType, water_velocity))

)

template<>
inline void getStateVector(const AdcpState &state, Eigen::Matrix<AdcpState::scalar, AdcpState::DOF, 1>& state_vector)
{
    state_vector.block(MTK::getStartIdx(&AdcpState::orientation),0,MTK::getDof(&AdcpState::orientation),1) = MTK::SO3<AdcpState::scalar>::log(state.orientation);
    state_vector.block(MTK::getStartIdx(&AdcpState::velocity),0,MTK::getDof(&AdcpState::velocity),1) = state.velocity;
    state_vector.block(MTK::getStartIdx(&AdcpState::position),0,MTK::getDof(&AdcpState::position),1) = state.position;
    state_vector.block(MTK::getStartIdx(&AdcpState::bias_gyro),0,MTK::getDof(&AdcpState::bias_gyro),1) = state.bias_gyro;
    state_vector.block(MTK::getStartIdx(&AdcpState::bias_acc),0,MTK::getDof(&AdcpState::bias_acc),1) = state.bias_acc;
    state_vector.block(MTK::getStartIdx(&AdcpState::bias_adcp),0,MTK::getDof(&AdcpState::bias_adcp),1) = state.bias_adcp;
    state_vector.block(MTK::getStartIdx(&AdcpState::water_velocity),0,MTK::getDof(&AdcpState::water_velocity),1) = state.water_velocity;
}

}

#endif
