#include "Mesh.hpp"
#include <Eigen/Core>
#include <memory>
#include "common.hpp"


#ifndef _MESH_CONVERT_HPP
#define _MESH_CONVERT_HPP

namespace PolyMesh2D
{
    namespace MeshTransform
    {
        std::unique_ptr<CurvedMesh::Mesh> mesh_convert(const StraightMesh::Mesh *s_mesh);
    }
}

#endif