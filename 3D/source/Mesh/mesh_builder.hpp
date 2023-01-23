#include "Mesh.hpp"
#include "mesh_reader.hpp"
#include <algorithm>
#include <memory> // for std::unique_pointer

#include <stdlib.h> /* exit, EXIT_FAILURE */

#ifndef MESH_BUILDER_HPP
#define MESH_BUILDER_HPP

namespace PolyMesh3D
{
    namespace StraightMesh
    {

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
            MeshBuilder(const std::string mesh_file);

            /**
             *  Build mesh
             */
            std::unique_ptr<Mesh> build_the_mesh();

        private:
            void build_boundary(Mesh *mesh);
            const std::string _mesh_file;
        };

    }
} // end namespace HArDCore3D
#endif
