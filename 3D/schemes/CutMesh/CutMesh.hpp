#include <Mesh.hpp>
#include <MeshBuilder2D.hpp>
#include <MeshReaderTyp2.hpp>
#include <function.hpp>
#include <common.hpp>
#include <memory>

#ifndef _CUTMESH_HPP
#define _CUTMESH_HPP

namespace PolyMesh2D
{
    class MeshCutter2
    {
    public:
        // MeshCutter2(StraightMesh::Mesh *s_mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param);
        MeshCutter2(CurvedMesh::Mesh *s_mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param);
        std::unique_ptr<CurvedMesh::Mesh> cut_mesh();

    private:
        void make_internal_cells(CurvedMesh::Mesh * c_mesh);
        void make_curved_cells(CurvedMesh::Mesh * c_mesh);
        void make_isotropic(CurvedMesh::Mesh * c_mesh);
        void make_boundary(CurvedMesh::Mesh * c_mesh);

        void make_convex(CurvedMesh::Mesh *c_mesh);

        // StraightMesh::Mesh *m_s_mesh;
        CurvedMesh::Mesh *m_s_mesh;
        Functional::ScalarFunction2D m_level_set;
        Functional::Curve m_param;

        bool no_external;
        bool make_straight;
    };

}

#endif