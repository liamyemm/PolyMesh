#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightFace.hpp"
#include "StraightCell.hpp"

namespace PolyMesh3D
{
    namespace StraightMesh
    {
        Edge::Edge(size_t index, std::array<Vertex *, 2> vertices)
            : Polytope::Polytope(index), _vertices(vertices)
        {
            _center_mass = 0.5 * (_vertices[0]->coords() + _vertices[1]->coords());
            _measure = (_vertices[1]->coords() - _vertices[0]->coords()).norm();
            _diameter = _measure;

            _tangent = (_vertices[1]->coords() - _vertices[0]->coords()) / _measure;
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

        std::array<VectorRd, 2> Edge::coords() const
        {
            return {_vertices[0]->coords(), _vertices[1]->coords()};
        }

        bool Edge::test() const
        {
            bool valid = true;
            if(_vertices.size() != 2)
            {
                std::cerr << "**** Edge " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has " << _vertices.size() << " vertices.\n";
                valid = false;
            }
            if(_measure < 1E-14)
            {
                std::cerr << "**** Edge " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has trivial measure: " << _measure << "\n";
                valid = false;
            }
            return valid;
        }
    }
}
