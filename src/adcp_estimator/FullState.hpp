#ifndef _FULL_STATE_HPP_
#define _FULL_STATE_HPP_

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

MTK_BUILD_MANIFOLD(FullState,
   ((RotationType, orientation))
   ((VelocityType, velocity)) // navigation/target frame velocity
   ((TranslationType, position)) // navigation/target frame velocity
   ((BiasType, bias_gyro))
   ((BiasType, bias_acc))
)

template<>
inline void getStateVector(const FullState &state, Eigen::Matrix<FullState::scalar, FullState::DOF, 1>& state_vector)
{
    state_vector.block(MTK::getStartIdx(&FullState::orientation),0,MTK::getDof(&FullState::orientation),1) = MTK::SO3<FullState::scalar>::log(state.orientation);
    state_vector.block(MTK::getStartIdx(&FullState::velocity),0,MTK::getDof(&FullState::velocity),1) = state.velocity;
    state_vector.block(MTK::getStartIdx(&FullState::position),0,MTK::getDof(&FullState::position),1) = state.position;
    state_vector.block(MTK::getStartIdx(&FullState::bias_gyro),0,MTK::getDof(&FullState::bias_gyro),1) = state.bias_gyro;
    state_vector.block(MTK::getStartIdx(&FullState::bias_acc),0,MTK::getDof(&FullState::bias_acc),1) = state.bias_acc;
}

}

#endif