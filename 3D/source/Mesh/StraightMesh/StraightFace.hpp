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

        class Face : public Polytope
        {
        public:
            Face(size_t index, std::vector<Edge *> edges);

            int edge_orientation(const size_t edge_index) const; ///< Return the orientation of the Edge located at edge_index
            
            VectorRd normal() const;      ///< Return the normal of the Face
            VectorRd edge_normal() const; ///< Return the normal of the Edge

            bool test() const override; ///< Return a boolean determining if the geometries of the Face are valid

        private:
            VectorRd _normal;

            std::vector<int> _edge_orientations;
            std::vector<VectorRd> _outer_normals;

            void set_edge_orientations();
        };

        //@}
    }
}

#endif