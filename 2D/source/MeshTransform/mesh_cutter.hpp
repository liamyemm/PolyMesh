#include "Mesh.hpp"

#include "function.hpp"
#include "common.hpp"

#include <memory>

#ifndef _MESH_CUTTER_HPP
#define _MESH_CUTTER_HPP

namespace PolyMesh2D
{
    // namespace MeshTransform
    // {
        class MeshCutter
        {
        public:
            MeshCutter(CurvedMesh::Mesh *mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param, bool keep_internal, bool keep_external);
            void cut_mesh();

        private:
            void make_uncut_cells(CurvedMesh::Mesh *c_mesh);
            void make_cut_cells(CurvedMesh::Mesh *c_mesh);
            void make_boundary(CurvedMesh::Mesh *c_mesh);

            CurvedMesh::Mesh *m_mesh;
            Functional::ScalarFunction2D m_level_set;
            Functional::Curve m_param;

            bool m_keep_internal;
            bool m_keep_external;
        };
    // }
}

#endif