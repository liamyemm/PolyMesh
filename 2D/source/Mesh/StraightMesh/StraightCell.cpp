#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightCell.hpp"

namespace PolyMesh2D
{
    namespace StraightMesh
    {
        double signed_area(VectorRd A, VectorRd B, VectorRd C) // returns signed area of trangle ABC. Positive if anticlock, negative otherwise
        {
            return (A(0) * (B(1) - C(1)) + B(0) * (C(1) - A(1)) + C(0) * (A(1) - B(1))) / 2.0;
        }

        VectorRd triangle_center_mass(Triangle triangle)
        {
            return (triangle[0] + triangle[1] + triangle[2]) / 3.0;
        }

        double triangle_measure(Triangle triangle)
        {
            VectorRd A = triangle[0];
            VectorRd B = triangle[1];
            VectorRd C = triangle[2];
            return std::abs((A(0) * (B(1) - C(1)) + B(0) * (C(1) - A(1)) + C(0) * (A(1) - B(1))) / 2.0);
        }

        Cell::Cell(size_t index, std::vector<Triangle> triangles) : Polytope::Polytope(index), _triangles(triangles)
        {
            std::vector<VectorRd> vertex_coords;

            for (auto &triangle : _triangles)
            {
                double measure_of_triangle = triangle_measure(triangle);
                _measure += measure_of_triangle;
                _center_mass += measure_of_triangle * triangle_center_mass(triangle);
                for (auto &coord : triangle)
                {
                    if (std::find(vertex_coords.begin(), vertex_coords.end(), coord) == vertex_coords.end())
                    {
                        vertex_coords.push_back(coord); // if coord not already in vertex_coords, add it
                    }
                }
            }
            _center_mass /= _measure;

            for (auto it = vertex_coords.begin(); it != vertex_coords.end(); ++it)
            {
                for (auto jt = it; jt != vertex_coords.end(); ++jt)
                {
                    _diameter = std::max(_diameter, (*it - *jt).norm());
                }
            }
            set_edge_orientations();
        }

        void Cell::set_edge_orientations() // not very efficient - probably room for improvement
        {
            for (size_t iF = 0; iF < _edges.size(); ++iF)
            {
                VectorRd normal = _edges[iF]->normal();
                std::array<VectorRd, 2> edge_coords = _edges[iF]->coords();

                VectorRd center;

                // find a cell simplex attached to the edge
                bool found = false;
                for (auto &cell_simplex : this->triangulation())
                {
                    double count = 0;
                    for (size_t i = 0; (i < cell_simplex.size()) && count < 2; ++i)
                    {
                        if (std::find(edge_coords.begin(), edge_coords.end(), cell_simplex[i]) == edge_coords.end()) // requires numerical precision
                        {
                            ++count;
                        }
                    }
                    if (count == 1) // only don't share one coordinate
                    {
                        center = triangle_center_mass(cell_simplex);
                        found = true;
                        break;
                    }
                }
                assert(found);

                _edge_orientations.push_back(Math::sgn((_edges[iF]->center_mass() - center).dot(normal)));
                assert(_edge_orientations[iF] != 0);
            }
        }

        int Cell::edge_orientation(const size_t edge_index) const
        {
            assert(edge_index < _edges.size());
            return (_edge_orientations[edge_index]);
        }

        VectorRd Cell::edge_normal(const size_t edge_index) const
        {
            assert(edge_index < _edges.size());
            return this->edge_orientation(edge_index) * (_edges[edge_index]->normal());
        }

        void Cell::plot_triangulation(std::ofstream *out) const
        {
            for (auto &triangle : _triangles)
            {
                for (size_t i = 0; i < 3; ++i)
                {
                    size_t i_next = (i + 1) % 3;
                    *out << triangle[i](0) << " " << triangle[i](1) << std::endl;
                    *out << triangle[i_next](0) << " " << triangle[i_next](1) << std::endl;
                    *out << std::endl;
                }
            }
        }
    }
}
