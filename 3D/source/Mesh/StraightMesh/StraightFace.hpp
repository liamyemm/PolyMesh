#include "StraightPolytope.hpp"
#include "function.hpp"

#ifndef _STRAIGHTFACE_HPP
#define _STRAIGHTFACE_HPP

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
            Edge(size_t index, std::array<VectorRd, 2> coords);

            int vertex_orientation(const size_t vertex_index) const; ///< Return the orientation of the Vertex located at vertex_index (1 if it is in the direction of the tangent, -1 otherwise)

            VectorRd tangent() const;                  ///< Return the tangent of the Edge
            VectorRd normal() const;                   ///< Return the normal of the Edge

            std::array<VectorRd, 2> coords() const; ///< Return an array of coordinates representing the coordinates of the vertices in the edge.

        private:
            std::array<VectorRd, 2> _coords;
            VectorRd _tangent;  
        };

        //@}
    }
}

#endif