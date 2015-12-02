#ifndef _POSE_WITH_VELOCITY_EXTRA_STATES_HPP_
#define _POSE_WITH_VELOCITY_EXTRA_STATES_HPP_

/** MTK library **/
#include <mtk/src/SubManifold.hpp>
#include <mtk/types/SOn.hpp>
#include <mtk/types/vect.hpp>
#include <mtk/startIdx.hpp>
#include <Eigen/Core>

namespace pose_estimation
{

template<class PoseWithVelocityType, class ExtraStatesType>
struct PoseWithVelocityExtraStates
{
    enum {PoseWithVelocityDOF = PoseWithVelocityType::DOF, ExtraStatesDOF = ExtraStatesType::DOF};
    enum {PoseWithVelocityIdx = 0, ExtraStatesIdx = PoseWithVelocityDOF};
	
public:
    typedef PoseWithVelocityType pose_with_velocity_type;
    typedef ExtraStatesType extrastates_type;
    typedef typename pose_with_velocity_type::scalar scalar;
    typedef PoseWithVelocityExtraStates self;
    enum {DOF = PoseWithVelocityDOF + ExtraStatesDOF};
    typedef Eigen::Matrix<unsigned, PoseWithVelocityDOF, 1> mask_type;

    // State types
    MTK::SubManifold<PoseWithVelocityType, PoseWithVelocityIdx> pose_with_velocity;
    MTK::SubManifold<ExtraStatesType, ExtraStatesIdx> extrastates;
    
    
    // Construct from pose and velocities
    PoseWithVelocityExtraStates(const pose_with_velocity_type &pose_with_velocity = pose_with_velocity_type(), const extrastates_type& extrastates = extrastates_type())
		    : pose_with_velocity(pose_with_velocity), extrastates(extrastates) {}
    
    // Copy constructor
    template<class PoseWithVelocityType2, class ExtraStatesType2>
    PoseWithVelocityExtraStates(const PoseWithVelocityExtraStates<PoseWithVelocityType2, ExtraStatesType2> &state)
		    : pose_with_velocity(state.pose_with_velocity), extrastates(state.extrastates) {}
    
    void boxplus(MTK::vectview<const scalar, DOF> state_vec, scalar _scale = 1)
    {
        pose_with_velocity.boxplus(MTK::subvector(state_vec, &PoseWithVelocityExtraStates::pose_with_velocity), _scale);
        extrastates.boxplus(MTK::subvector(state_vec, &PoseWithVelocityExtraStates::extrastates), _scale);
    }
    void boxminus(MTK::vectview<scalar, DOF> state_ret, const PoseWithVelocityExtraStates &other) const
    {
        pose_with_velocity.boxminus(MTK::subvector(state_ret, &PoseWithVelocityExtraStates::pose_with_velocity), other.pose_with_velocity);
        extrastates.boxminus(MTK::subvector(state_ret, &PoseWithVelocityExtraStates::extrastates), other.extrastates);
    }
    
    friend std::ostream& operator<<(std::ostream &os, const PoseWithVelocityExtraStates<PoseWithVelocityType, ExtraStatesType> &state)
    {
	return os << state.pose_with_velocity << " " << state.extrastates;
    }
    friend std::istream& operator>>(std::istream &is, PoseWithVelocityExtraStates<PoseWithVelocityType, ExtraStatesType> &state)
    {
	return is >> state.pose_with_velocity >> state.extrastates;
    }

    void applyVelocity(double delta_time)
    {
        pose_with_velocity.applyVelocity(delta_time);
    }

    void applyVelocity(const Eigen::Vector3d& acc, double delta_time)
    {
        pose_with_velocity.applyVelocity(acc, delta_time);
    }
    
    Eigen::Matrix<scalar, PoseWithVelocityDOF, 1> getStateVector() const
    {
        return pose_with_velocity.getStateVector();
    }

    Eigen::Matrix<scalar, -1, 1> getSubStateVector(const mask_type& mask) const
    {
        return pose_with_velocity.getSubStateVector(mask);
    }
};

}

#endif