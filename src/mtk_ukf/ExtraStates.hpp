#ifndef _EXTRA_STATES_HPP_
#define _EXTRA_STATES_HPP_

/** MTK library **/
#include <mtk/src/SubManifold.hpp>
#include <mtk/types/vect.hpp>
#include <mtk/startIdx.hpp>
#include <Eigen/Core>

namespace pose_estimation
{

template<int dof>
struct ExtraStates
{
    enum {ExtraStatesDOF = MTK::vect<dof>::DOF};
    enum {ExtraStatesIdx = 0};
public:
    typedef MTK::vect<dof> extrastates_type;
    typedef typename extrastates_type::scalar scalar;
    typedef ExtraStates self;
    enum {DOF = ExtraStatesDOF};
    typedef Eigen::Matrix<unsigned, DOF, 1> map_type;

    // State types
    MTK::SubManifold<MTK::vect<dof>, ExtraStatesIdx> extrastates;

    // Helper
    map_type extrastates_map;

    // Construct from pose and velocities
    ExtraStates(const extrastates_type& extrastates = extrastates_type()) : extrastates(extrastates), extrastates_map(map_type::Zero()) {}

    // Copy constructor
    template<int dof2>
    ExtraStates(const ExtraStates<dof2> &state) : extrastates(state.extrastates), extrastates_map(map_type::Zero()) {}

    void boxplus(MTK::vectview<const scalar, DOF> state_vec, scalar _scale = 1)
    {
        extrastates.boxplus(MTK::subvector(state_vec, &ExtraStates::extrastates), _scale);
    }
    void boxminus(MTK::vectview<scalar, DOF> state_ret, const ExtraStates &other) const
    {
        extrastates.boxminus(MTK::subvector(state_ret, &ExtraStates::extrastates), other.extrastates);
    }

    friend std::ostream& operator<<(std::ostream &os, const ExtraStates<dof> &state)
    {
        return os << state.extrastates;
    }
    friend std::istream& operator>>(std::istream &is, ExtraStates<dof> &state)
    {
        return is >> state.extrastates;
    }

    void setExtraStatesMap(const map_type& map)
    {
        extrastates_map = map;

        // make sure the indeces in the map are ordered
        std::sort(extrastates_map.data(), extrastates_map.data()+extrastates_map.size());
    }

    map_type getExtraStatesMap()
    {
        return extrastates_map;
    }

    template<int mask_dim>
    Eigen::Matrix<scalar, -1, 1> getExtraStates(const Eigen::Matrix<unsigned, mask_dim, 1>& mask) const
    {
        Eigen::Matrix<scalar, -1, 1> extrastates_state;
        extrastates_state.resize(mask.count());
        extrastates_state.setZero();

        // populate extra states sub state
        unsigned sub_state_offset = 0;
        for(unsigned i = 0; i < mask_dim; i++)
        {
            if(mask[i] > 0)
            {
                for(unsigned j = 0; j < DOF; j++)
                {
                    if(extrastates_map[j] == i)
                    {
                        extrastates_state[sub_state_offset] = extrastates[j];
                    }
                }
                sub_state_offset++;
            }
        }

        return extrastates_state;
    }
};

}

#endif