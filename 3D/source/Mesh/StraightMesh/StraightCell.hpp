#include "StraightPolytope.hpp"

#ifndef _STRAIGHTCELL_HPP
#define _STRAIGHTCELL_HPP

namespace PolyMesh3D
{
    namespace StraightMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Cell : public Polytope
        {
        public:
            Cell(size_t index, std::vector<Face *> faces); ///< Constructor for a Cell. Calls the base class constructor, sets _triangles to @arg triangles and computes the cell data (center_mass, diameter, measure, edge orientations)

            int face_orientation(const size_t face_index) const; ///< Return the orientation of the Face located at @arg face_index.
            VectorRd face_normal(const size_t face_index) const; ///< Return the outer normal of the Face located at @arg face_index.

            bool test() const override; ///< Return a boolean determining if the geometries of the Cell are valid

        private:
            std::vector<int> _face_orientations;
            void set_face_orientations();
        };

        //@}
    }
}

#endif