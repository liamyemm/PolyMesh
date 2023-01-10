#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightFace.hpp"
#include "StraightCell.hpp"

namespace PolyMesh3D
{
    namespace StraightMesh
    {
        Vertex::Vertex(size_t index, VectorRd coords) : Polytope::Polytope(index), _coords(coords)
        {
            _measure = 1.0;
            _center_mass = _coords;
        }

        VectorRd Vertex::coords() const
        {
            return _coords;
        }

        bool Vertex::test() const
        {
            bool valid = true;
            if(_vertices.size() != _edges.size())
            {
                std::cerr << "**** Vertex " << _index <<" located at (" << _coords(0) << ", " << _coords(1) << ", " << _coords(2) << ") has " << _vertices.size() << " vertices and " << _edges.size() << " edges.\n";
                valid = false;
            }
            return valid;
        }
    }
}