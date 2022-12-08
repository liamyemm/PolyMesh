#include <iostream>

#include "Vertex.hpp"
#include "Edge.hpp"
#include "Cell.hpp"

#include "GaussLegendre.hpp"
#include "function.hpp"

namespace PolyMesh2D
{
    namespace CurvedMesh
    {
        Cell::Cell(size_t index, std::vector<Edge *> edges) : Polytope::Polytope(index)
        {
            _edges = edges;
            set_edge_orientations();

            assert(_edge_orientations.size() == _edges.size());

            _is_straight = true;
            for (auto &e : _edges)
            {
                _is_straight = (_is_straight && e->is_straight());
            }

            unsigned partitions = (_is_straight ? 1 : 10);

            for (size_t i = 0; i < _edges.size(); ++i)
            {
                double deltat_i = (_edges[i]->parameterisation().tmax - _edges[i]->parameterisation().tmin) / double(partitions);
                for (unsigned part_i = 0; part_i < partitions; ++part_i)
                {
                    double t_val_i = _edges[i]->parameterisation().tmin + part_i * deltat_i;
                    VectorRd point_i = _edges[i]->parameterisation().value(t_val_i);
                    for (size_t j = 0; j < _edges.size(); ++j)
                    {
                        double deltat_j = (_edges[j]->parameterisation().tmax - _edges[j]->parameterisation().tmin) / double(partitions);
                        for (unsigned part_j = 0; part_j < partitions; ++part_j)
                        {
                            double t_val_j = _edges[j]->parameterisation().tmin + part_j * deltat_j;
                            VectorRd point_j = _edges[j]->parameterisation().value(t_val_j);
                            _diameter = std::max(_diameter, (point_i - point_j).norm());
                        }
                    }
                }
            }

            unsigned quad_degree = (_is_straight ? 2 : 20); // if straight, max degree of integration is that of x(x.n) which is 2.

            Quadrature::GaussLegendre1D quad(quad_degree);

            // compute _measure and _center_mass via integration on boundary
            for (size_t iE = 0; iE < _edges.size(); ++iE)
            {
                Functional::Curve edge_param = _edges[iE]->parameterisation();

                double a = edge_param.tmin;
                double b = edge_param.tmax;

                double edge_len = b - a;
                assert(edge_len > 0.0);

                for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
                {
                    double point = edge_len * quad.point(iqn) + a;
                    double tmp = edge_len * quad.weight(iqn) * edge_param.value(point).dot(this->edge_normal(iE, point)) * edge_param.derivative(point).norm();
                    _measure += tmp;
                    _center_mass += edge_param.value(point) * tmp;
                }
            }
            _measure = _measure / 2.0;
            _center_mass = _center_mass / (3.0 * _measure);
        }

        void Cell::set_edge_orientations()
        {
            _edge_orientations.clear(); // reset to empty vector
            VectorRd temp_center = Eigen::Vector2d::Zero();
            for (auto &e : _edges)
            {
                temp_center += e->parameterisation().value(0.5 * (e->parameterisation().tmin + e->parameterisation().tmax));
            }
            temp_center = temp_center / _edges.size();

            for (size_t iE = 0; iE < _edges.size(); ++iE)
            {
                Edge *edge = _edges[iE];
                VectorRd normal_at_edge_center = edge->parameterisation().normal(0.5 * (edge->parameterisation().tmin + edge->parameterisation().tmax));
                // VectorRd value_at_edge_center = edge->parameterisation().value(0.5 * (edge->parameterisation().tmin + edge->parameterisation().tmax));
                _edge_orientations.push_back(Math::sgn(normal_at_edge_center.dot(edge->center_mass() - temp_center)));
            }
        }

        int Cell::edge_orientation(const size_t edge_index) const
        {
            assert(edge_index < _edges.size());
            return (_edge_orientations[edge_index]);
        }

        VectorRd Cell::edge_normal(const size_t edge_index, double t) const
        {
            assert(edge_index < _edges.size());
            return this->edge_orientation(edge_index) * (_edges[edge_index]->normal(t));
        }

        bool Cell::is_straight() const
        {
            return _is_straight;
        }

        bool Cell::test() const
        {
            bool valid = true;
            if (_vertices.size() != _edges.size())
            {
                std::cout << "Error! Cell " << _index << " has " << _vertices.size() << " vertices and " << _edges.size() << " edges.\n";
                valid = false;
            }

            double integral = 0.0;
            Eigen::Vector2d const_vec(1.0, 1.0);

            unsigned quad_degree = (_is_straight ? 0 : 20); // if straight, max degree of integration is 0.

            Quadrature::GaussLegendre1D quad(quad_degree);

            // compute _measure and _center_mass via integration on boundary
            for (size_t iE = 0; iE < _edges.size(); ++iE)
            {
                Functional::Curve edge_param = _edges[iE]->parameterisation();

                double a = edge_param.tmin;
                double b = edge_param.tmax;

                double edge_len = b - a;
                assert(edge_len > 0.0);

                for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
                {
                    double point = edge_len * quad.point(iqn) + a;
                    integral += edge_len * quad.weight(iqn) * const_vec.dot(this->edge_normal(iE, point)) * edge_param.derivative(point).norm();
                }
            }

            if (std::abs(integral) > 1E-12)
            {
                std::cout << "Error! Cell " << _index << " has boundary integral evaluated to " << integral << ".\n";
                valid = false;
            }
            return valid;
        }
    }
}

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
