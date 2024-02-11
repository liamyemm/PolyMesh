#include "affine_shift.hpp"
#include "MeshReaderTyp2.hpp"

#include <iomanip>

namespace PolyMesh2D
{
    void affine_shift(const std::string &input_file_path, const std::string &output_file_path, const Eigen::Matrix2d &A, const Eigen::Vector2d &b)
    {
        PolyMesh2D::StraightMesh::MeshReaderTyp2 mesh_reader(input_file_path);

        std::vector<std::array<double, 2>> vertices;
        std::vector<std::vector<size_t>> cells;

        mesh_reader.read_mesh(vertices, cells);

        std::ofstream out(output_file_path);

        out << "Vertices"
            << "\n";
        out << vertices.size() << "\n";
        for (size_t iV = 0; iV < vertices.size(); iV++)
        {
            Eigen::Vector2d v;
            v << vertices[iV][0], vertices[iV][1];
            v = A * v + b;
            out << std::setprecision(17) << v(0);
            out << " ";
            out << std::setprecision(17) << v(1);
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
    }
}
