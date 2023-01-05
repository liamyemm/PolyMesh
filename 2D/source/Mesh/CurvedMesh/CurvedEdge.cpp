#include <iostream>

#include "CurvedVertex.hpp"
#include "CurvedEdge.hpp"
#include "CurvedCell.hpp"

#include "GaussLegendre.hpp"

namespace PolyMesh2D
{
    namespace CurvedMesh
    {

        Edge::Edge(size_t index, Functional::Curve &parameterisation)
            : Polytope::Polytope(index), _parameterisation(parameterisation)
        {
            // compute _measure and _center_mass via integration
            Quadrature::GaussLegendre1D quad(20);

            double a = _parameterisation.tmin;
            double b = _parameterisation.tmax;

            double len = b - a;
            assert(len > 0.0);

            for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
            {
                double point = len * quad.point(iqn) + a;
                double tmp = len * quad.weight(iqn) * _parameterisation.derivative(point).norm();
                _measure += tmp;
                _center_mass += _parameterisation.value(point) * tmp; // is this really what we want? It is not necessarily on the curve
            }
            _center_mass /= _measure;
            _diameter = _measure; // not actually true, but will do for now.

            _is_straight = false; // initialise to false and set true if edge is straight
        }

        int Edge::vertex_orientation(const size_t vertex_index) const
        {
            assert(vertex_index < _vertices.size());
            VectorRd tangent_at_center = this->tangent(0.5 * (_parameterisation.tmin + _parameterisation.tmax));
            return Math::sgn(tangent_at_center.dot(_vertices[vertex_index]->coords() - _center_mass));
        }

        VectorRd Edge::tangent(double t) const
        {
            return _parameterisation.tangent(t);
        }

        VectorRd Edge::normal(double t) const
        {
            return _parameterisation.normal(t);
        }

        const Functional::Curve &Edge::parameterisation() const
        {
            return _parameterisation;
        }

        std::array<VectorRd, 2> Edge::coords() const
        {
            return {_parameterisation.value(_parameterisation.tmin), _parameterisation.value(_parameterisation.tmax)};
        }

        void Edge::plot(std::ofstream *out, int partitions) const
        {
            if (_is_straight)
            {
                partitions = 1;
            }
            double deltat = (_parameterisation.tmax - _parameterisation.tmin) / double(partitions);
            for (int i = 0; i < partitions - 1; ++i)
            {
                double t0 = _parameterisation.tmin + i * deltat;
                double t1 = _parameterisation.tmin + (i + 1) * deltat;

                *out << _parameterisation.value(t0)(0) << " " << _parameterisation.value(t0)(1) << std::endl;
                *out << _parameterisation.value(t1)(0) << " " << _parameterisation.value(t1)(1) << std::endl;
                *out << std::endl;
            }
            double t0 = _parameterisation.tmax - deltat;
            double t1 = _parameterisation.tmax;

            *out << _parameterisation.value(t0)(0) << " " << _parameterisation.value(t0)(1) << std::endl;
            *out << _parameterisation.value(t1)(0) << " " << _parameterisation.value(t1)(1) << std::endl;
            *out << std::endl;
        }

        bool Edge::is_straight() const
        {
            return _is_straight;
        }

        void Edge::set_straight()
        {
            _is_straight = true;
        }

        bool Edge::test() const
        {
            bool valid = true;
            if (_vertices.size() != 2)
            {
                std::cout << "Error! Edge " << _index << " has " << _vertices.size() << " vertices.\n";
                valid = false;
            }
            if (is_boundary())
            {
                if (_cells.size() != 1)
                {
                    std::cout << "Error! Edge " << _index << " is a boundary edge and has " << _cells.size() << " cells.\n";
                    valid = false;
                }
            }
            else
            {
                if (_cells.size() != 2)
                {
                    std::cout << "Error! Edge " << _index << " is not a boundary edge and has " << _cells.size() << " cells.\n";
                    valid = false;
                }
            }

            Eigen::Vector2d start_point(_parameterisation.value(_parameterisation.tmin));
            Eigen::Vector2d end_point(_parameterisation.value(_parameterisation.tmax));

            Eigen::Vector2d v0(_vertices[0]->coords());
            Eigen::Vector2d v1(_vertices[1]->coords());

            if (((start_point - v0).norm() > 1E-12) || ((end_point - v1).norm() > 1E-12))
            {
                if (((start_point - v1).norm() > 1E-12) || ((end_point - v0).norm() > 1E-12))
                {
                    std::cout << "Error! Edge " << _index << " has inconsistent vertices.\n";
                    valid = false;
                }
            }

            return valid;
        }
    }
}