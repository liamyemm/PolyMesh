#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightCell.hpp"

namespace PolyMesh2D
{
    namespace StraightMesh
    {
        Vertex::Vertex(size_t index, VectorRd coords)
             : Polytope::Polytope(index), _coords(coords)
        {
            _measure = 1.0;
            _center_mass = _coords;
        }

        VectorRd Vertex::coords() const
        {
            return _coords;
        }
    }
}