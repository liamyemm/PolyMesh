#include "CurvedPolytope.hpp"

#ifndef _CURVEDCELL_HPP
#define _CURVEDCELL_HPP

namespace PolyMesh2D
{
    namespace CurvedMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Cell : public Polytope
        {
        public:
            Cell(size_t index, std::vector<Edge *> edges); ///< Constructor for a Cell. Calls the base class constructor, sets _edges to @arg edges and computes the cell data (center_mass, diameter, measure, edge orientations)

            int edge_orientation(const size_t edge_index) const;           ///< Return the orientation of the Edge located at @arg edge_index.
            VectorRd edge_normal(const size_t edge_index, double t) const; ///< Return the outer normal of the Edge located at @arg edge_index evaluated at the parameter @t.

            bool is_straight() const; ///< Return a boolean stating whether or not all the edges of the Cell are straight.

            bool test() const;

            void set_edge_orientations();

        private:
            bool _is_straight;
            std::vector<int> _edge_orientations;
        };

        //@}
    }
}

#endif