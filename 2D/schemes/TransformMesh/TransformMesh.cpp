#include <Mesh.hpp>
#include <TransformMesh.hpp>
#include <iomanip>
#include <MeshReaderTyp2.hpp>

int main(int arc, char **argv)
{
    PolyMesh2D::StraightMesh::MeshReaderTyp2 mesh_reader("../../../typ2_meshes/mesh2_6.typ2");

    std::vector<std::array<double, 2>> vertices;
    std::vector<std::vector<size_t>> cells;

    mesh_reader.read_mesh(vertices, cells);

    // std::vector<Eigen::Vector2d> new_verts;

    // for (size_t iV = 0; iV < vertices.size(); ++iV)
    // {
    //     double v0 = 2.0 * vertices[iV][0] - 1.0;
    //     double v1 = 2.0 * vertices[iV][1] - 1.0;
    //     new_verts.push_back(Eigen::Vector2d(v0, v1));
    // }

    std::ofstream out("../../../typ2_meshes/mesh2_6_transformed.typ2");

    out << "Vertices"
        << "\n";
    out << vertices.size() << "\n";
    for (size_t iV = 0; iV < vertices.size(); iV++)
    {
        out << std::setprecision(17) << (2.0 * vertices[iV][0] - 1.0);
        out << " ";
        out << std::setprecision(17) << (2.0 * vertices[iV][1] - 1.0);
        out << "\n";
    }

    out << "cells"
        << "\n";
    out << cells.size() << "\n";
    for (std::size_t cell_i = 0; cell_i < cells.size(); cell_i++)
    {
        out << std::setw(6) << cells[cell_i].size();
        for (std::size_t node_j = 0; node_j < cells[cell_i].size(); node_j++)
        {
            out << std::setw(6) << cells[cell_i][node_j];
        }
        out << "\n";
    }

    out.close();

    return 0;
}
