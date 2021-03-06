#ifndef _POSE_WITH_VELOCITY_HPP_
#define _POSE_WITH_VELOCITY_HPP_

/** MTK library **/
#include <mtk/src/SubManifold.hpp>
#include <mtk/types/SOn.hpp>
#include <mtk/types/vect.hpp>
#include <mtk/build_manifold.hpp>
#include <mtk/startIdx.hpp>
#include <boost/concept_check.hpp>
#include <pose_estimation/ManifoldHelper.hpp>

namespace pose_estimation
{

typedef MTK::SO3<double> RotationType;
typedef RotationType::vect_type TranslationType;
typedef RotationType::vect_type VelocityType;

MTK_BUILD_MANIFOLD(PoseWithVelocity,
   ((TranslationType, position))
   ((RotationType, orientation))
   ((VelocityType, velocity))
   ((VelocityType, angular_velocity))
)

template<>
inline void getStateVector(const PoseWithVelocity &state, Eigen::Matrix<PoseWithVelocity::scalar, PoseWithVelocity::DOF, 1>& state_vector)
{
    state_vector.block(MTK::getStartIdx(&PoseWithVelocity::position),0,MTK::getDof(&PoseWithVelocity::position),1) = state.position;
    state_vector.block(MTK::getStartIdx(&PoseWithVelocity::orientation),0,MTK::getDof(&PoseWithVelocity::orientation),1) = MTK::SO3<PoseWithVelocity::scalar>::log(state.orientation);
    state_vector.block(MTK::getStartIdx(&PoseWithVelocity::velocity),0,MTK::getDof(&PoseWithVelocity::velocity),1) = state.velocity;
    state_vector.block(MTK::getStartIdx(&PoseWithVelocity::angular_velocity),0,MTK::getDof(&PoseWithVelocity::angular_velocity),1) = state.angular_velocity;
}

}

#endif