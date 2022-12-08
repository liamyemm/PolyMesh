#include <Mesh.hpp>
#include <MeshBuilder2D.hpp>
#include <MeshReaderTyp2.hpp>
#include <function.hpp>
#include <common.hpp>
#include <memory>

#ifndef _INTERSECTMESH_HPP
#define _INTERSECTMESH_HPP

namespace PolyMesh2D
{
    class MeshCutter
    {
    public:
        MeshCutter(StraightMesh::Mesh *s_mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param);
        std::unique_ptr<CurvedMesh::Mesh> cut_mesh();

    private:
        void make_internal_cells(CurvedMesh::Mesh * c_mesh);
        void make_convex(CurvedMesh::Mesh * c_mesh);
        void make_curved_cells(CurvedMesh::Mesh * c_mesh);
        void make_isotropic(CurvedMesh::Mesh * c_mesh);
        void make_homogeneous(CurvedMesh::Mesh * c_mesh);
        void make_boundary(CurvedMesh::Mesh * c_mesh);

        StraightMesh::Mesh *m_s_mesh;
        Functional::ScalarFunction2D m_level_set;
        Functional::Curve m_param;
    };

}

#endif