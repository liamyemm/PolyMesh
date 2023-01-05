#include <iostream>

#include "CurvedVertex.hpp"
#include "CurvedEdge.hpp"
#include "CurvedCell.hpp"

namespace PolyMesh2D
{
    namespace CurvedMesh
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


        bool Vertex::test() const
        {
            bool valid = true;
            if(_vertices.size() != _edges.size())
            {
                std::cout << "Error! Vertex " << _index << " has " << _vertices.size() << " vertices and " << _edges.size() << " edges.\n";
                valid = false;
            }
            if(!_is_boundary)
            {
                if(_vertices.size() != _cells.size())
                {
                    std::cout << "Error! Internal Vertex " << _index << " has " << _vertices.size() << " vertices and " << _cells.size() << " cells.\n";
                    valid = false;
                }
            }
            else
            {
                if(_vertices.size() != _cells.size() + 1)
                {
                    std::cout << "Error! Boundary Vertex " << _index << " has " << _vertices.size() << " vertices and " << _cells.size() << " cells.\n";
                    valid = false;
                }
            }
            return valid;
        }
    }
}