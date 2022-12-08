#include <Polytope.hpp>

#ifndef _VERTEX_HPP
#define _VERTEX_HPP

namespace PolyMesh2D
{
    namespace CurvedMesh
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

            bool test() const;

        private:
            VectorRd _coords;
        };

        //@}
    }

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

        private:
            VectorRd _coords;
        };

        //@}
    }
}

#endif