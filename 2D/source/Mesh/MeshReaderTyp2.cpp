#include "MeshReaderTyp2.hpp"

namespace PolyMesh2D
{
    namespace StraightMesh
    {
        MeshReaderTyp2::MeshReaderTyp2(std::string file_name) : _file_name(file_name) {}

        void MeshReaderTyp2::read_mesh(std::vector<std::array<double, 2>> &vertices, std::vector<std::vector<std::size_t>> &cells)
        {
            std::ifstream inFile;

            std::cout << "[MeshReaderTyp2] Reading mesh file " + _file_name + "\n";

            inFile.open(_file_name);
            if (!inFile)
            {
                std::cerr << "     Unable to open mesh file\n";
                exit(1);
            }

            // std::cout << "     Reading vertices...\n";

            //  Read mesh file
            std::string ignore_line;
            std::getline(inFile, ignore_line); // ignore the first line (should simply say "Vertices")

            std::size_t n_verts;
            inFile >> n_verts;

            std::size_t vert_count = 0;
            double xcoord, ycoord;
            std::array<double, 2> coord;

            while (vert_count < n_verts && inFile >> xcoord >> ycoord)
            {
                coord[0] = xcoord;
                coord[1] = ycoord;
                vertices.push_back(coord);
                vert_count++;
            }

            std::cout << "    Read " + std::to_string(vertices.size()) + "/" + std::to_string(n_verts) + " vertices\n";

            std::getline(inFile, ignore_line); // ignore the remainder of previous line
            std::getline(inFile, ignore_line); // ignore the following line (should simply say "cells")

            std::size_t n_cells;
            inFile >> n_cells;

            std::size_t cell_count = 0;
            std::size_t n_local_cell_nodes;

            while (cell_count < n_cells && inFile >> n_local_cell_nodes)
            {
                std::size_t node_count = 0;
                std::vector<std::size_t> cell_nodes(n_local_cell_nodes);
                while (node_count < n_local_cell_nodes && inFile >> cell_nodes[node_count])
                {
                    node_count++;
                }
                cells.push_back(cell_nodes);
                cell_count++;
            }

            std::cout << "    Read " + std::to_string(cells.size()) + "/" + std::to_string(n_cells) + " cells\n\n";

            inFile.close();
        }
    }
}
