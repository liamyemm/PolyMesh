
#include <memory> //std::unique_ptr
#include "MeshReaderTyp2.hpp"
#include "Polytope.hpp"
#include "Mesh.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"
#include "Cell.hpp"

#ifndef MESHBUILDER2D_HPP
#define MESHBUILDER2D_HPP

namespace PolyMesh2D
{
    namespace StraightMesh
    {

        /*!
         *  @addtogroup Mesh
         * @{
         */

        // ----------------------------------------------------------------------------
        //                            Class definition
        // ----------------------------------------------------------------------------

        /// The MeshBuilder class provides build tools to create a full mesh with all connectivities
        class MeshBuilder
        {
        public:
            /**
             * Constructor for MeshBuilder.
             */
            MeshBuilder() {}

            /**
             * Overloaded constructor for MeshBuilder so read_mesh can be called from build_the_mesh().
             */
            MeshBuilder(const std::string mesh_file);

            /**
             *  Build mesh
             */
            std::unique_ptr<Mesh> build_the_mesh();

        private:
            void build_boundary(Mesh *mesh); ///< identifies boundary cells and vertices, and compute lists of boundary cells, edges and vertices
            const std::string _mesh_file;
        };

        /*@}*/
    }
}
#endif /* MESHBUILDER2D_HPP */
