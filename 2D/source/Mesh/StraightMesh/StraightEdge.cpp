#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightCell.hpp"

#include "GaussLegendre.hpp"

namespace PolyMesh2D
{
    namespace StraightMesh
    {
        Edge::Edge(size_t index, std::array<VectorRd, 2> coords)
            : Polytope::Polytope(index), _coords(coords)
        {
            _center_mass = 0.5 * (_coords[0] + _coords[1]);
            _measure = (_coords[1] - _coords[0]).norm();
            _diameter = _measure;

            _tangent = (_coords[1] - _coords[0]) / _measure;
            _normal = VectorRd(-_tangent(1), _tangent(0));
        }

        int Edge::vertex_orientation(const size_t vertex_index) const
        {
            assert(vertex_index < _vertices.size());
            return Math::sgn(_tangent.dot(_vertices[vertex_index]->coords() - _center_mass));
        }

        VectorRd Edge::tangent() const
        {
            return _tangent;
        }

        VectorRd Edge::normal() const
        {
            return _normal;
        }

        std::array<VectorRd, 2> Edge::coords() const
        {
            return _coords;
        }

        void Edge::plot(std::ofstream *out) const
        {
            *out << _coords[0](0) << " " << _coords[0](1) << std::endl;
            *out << _coords[1](0) << " " << _coords[1](1) << std::endl;
            *out << std::endl;
        }
    }
}
