#include <iostream>

#include "CurvedVertex.hpp"
#include "CurvedEdge.hpp"
#include "CurvedCell.hpp"

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
                    double t_val_i_1 = _edges[i]->parameterisation().tmin + part_i * deltat_i;
                    double t_val_i_2 = _edges[i]->parameterisation().tmax - part_i * deltat_i;
                    VectorRd point_i_1 = _edges[i]->parameterisation().value(t_val_i_1);
                    VectorRd point_i_2 = _edges[i]->parameterisation().value(t_val_i_2);
                    for (size_t j = 0; j < _edges.size(); ++j)
                    {
                        double deltat_j = (_edges[j]->parameterisation().tmax - _edges[j]->parameterisation().tmin) / double(partitions);
                        for (unsigned part_j = 0; part_j < partitions; ++part_j)
                        {
                            double t_val_j_1 = _edges[j]->parameterisation().tmin + part_j * deltat_j;
                            double t_val_j_2 = _edges[j]->parameterisation().tmax - part_j * deltat_j;
                            VectorRd point_j_1 = _edges[j]->parameterisation().value(t_val_j_1);
                            VectorRd point_j_2 = _edges[j]->parameterisation().value(t_val_j_2);
                            _diameter = std::max(_diameter, (point_i_1 - point_j_1).norm());
                            _diameter = std::max(_diameter, (point_i_1 - point_j_2).norm());
                            _diameter = std::max(_diameter, (point_i_2 - point_j_1).norm());
                            _diameter = std::max(_diameter, (point_i_2 - point_j_2).norm());
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
