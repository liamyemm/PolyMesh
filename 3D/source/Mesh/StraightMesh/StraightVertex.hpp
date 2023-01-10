#include "StraightPolytope.hpp"

#ifndef _VERTEX_HPP
#define _VERTEX_HPP

namespace PolyMesh3D
{
    namespace StraightMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Vertex : public Polytope
        {
        public:
            Vertex(size_t, VectorRd);

            VectorRd coords() const; ///< Return the coordinates of a Vertex

            bool test() const override; ///< Return a boolean determining if the geometries of the Vertex are valid

        private:
            VectorRd _coords;
        };

        //@}
    }
}

#endif