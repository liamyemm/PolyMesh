#ifndef XVEMCORE_H
#define XVEMCORE_H

#include <gtest/gtest.h>
#include "../XVEMCore.h"

#include "MeshBuilder2D.hpp"
#include "MeshReaderTyp2.hpp"

#include "mesh_convert.hpp"

namespace PolyMesh2D
{
    class XVEMCoreTest : public ::testing::Test
    {
    protected:
        std::unique_ptr<PolyMesh2D::CurvedMesh::Mesh> curved_mesh;

        void SetUp() override
        {
            // Initialize the mesh object before each test case
            // Build mesh and reorder edges
            // PolyMesh2D::StraightMesh::MeshBuilder builder = PolyMesh2D::StraightMesh::MeshBuilder("/home/liam/github/codes/liamyemm/PolyMesh/2D/typ2_meshes/hexa1_2.typ2");
            PolyMesh2D::StraightMesh::MeshBuilder builder = PolyMesh2D::StraightMesh::MeshBuilder("/home/liam/github/codes/liamyemm/PolyMesh/2D/typ2_meshes/mesh2_1.typ2");
            std::unique_ptr<PolyMesh2D::StraightMesh::Mesh> straight_mesh;

            try
            {
                straight_mesh = builder.build_the_mesh();
            }
            catch (std::string msg)
            {
                std::cerr << msg;
            }

            curved_mesh = PolyMesh2D::MeshTransform::mesh_convert(straight_mesh.get());

            straight_mesh.reset(); // delete straight mesh
        }

        // void TearDown() override { /* Optional teardown */ }
    };
}

#endif // XVEMCORE_H