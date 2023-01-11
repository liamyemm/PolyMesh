#include "StraightPolytope.hpp"

#ifndef _STRAIGHTEDGE_HPP
#define _STRAIGHTEDGE_HPP

namespace PolyMesh2D
{
    namespace StraightMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Edge : public Polytope
        {
        public:
            Edge(size_t index, std::array<Vertex *, 2> vertices);

            int vertex_orientation(const size_t vertex_index) const; ///< Return the orientation of the Vertex located at vertex_index (1 if it is in the direction of the tangent, -1 otherwise)

            VectorRd tangent() const;                  ///< Return the tangent of the Edge
            std::array<VectorRd, 2> coords() const; ///< Return an array of coordinates representing the coordinates of the vertices in the edge.

            bool test() const override; ///< Return a boolean determining if the geometries of the Edge are valid

        private:
            VectorRd _tangent;  
        };

        //@}
    }
}

#endif