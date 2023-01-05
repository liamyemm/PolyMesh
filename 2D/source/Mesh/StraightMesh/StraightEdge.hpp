#include "StraightPolytope.hpp"
#include "function.hpp"

#ifndef _EDGE_HPP
#define _EDGE_HPP

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

            int vertex_orientation(const size_t vertex_index) const; ///< Return the orientation of the Vertex located at vertex_index

            VectorRd tangent() const;                  ///< Return the tangent of the Edge at the parameter @arg t
            VectorRd normal() const;                   ///< Return the normal of the Edge at the parameter @arg t

            std::array<VectorRd, 2> coords() const;

            void plot(std::ofstream *out) const;

        private:
            std::array<VectorRd, 2> _coords;

            VectorRd _normal;
            VectorRd _tangent;  
        };

        //@}
    }
}

#endif