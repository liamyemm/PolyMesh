#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightFace.hpp"
#include "StraightCell.hpp"

namespace PolyMesh3D
{
    namespace StraightMesh
    {
        Cell::Cell(size_t index, std::vector<Face *> faces) : Polytope::Polytope(index)
        {
            _faces = faces;
            std::vector<VectorRd> vertex_coords;
            for(auto & f : faces)
            {
                for(auto &v : f->get_vertices())
                {
                    if (std::find(vertex_coords.begin(), vertex_coords.end(), v->coords()) == vertex_coords.end())
                    {
                        vertex_coords.push_back(v->coords()); // if coord not already in vertex_coords, add it
                    }
                }
            }

            for (auto it = vertex_coords.begin(); it != vertex_coords.end(); ++it)
            {
                for (auto jt = it; jt != vertex_coords.end(); ++jt)
                {
                    _diameter = std::max(_diameter, (*it - *jt).norm());
                }
            }

            set_face_orientations();
            
            for (size_t iF = 0; iF < _faces.size(); ++iF)
            {
                double temp = _faces[iF]->measure() * _faces[iF]->center_mass().dot(this->face_normal(iF));
                _measure += temp;
                _center_mass += temp * _faces[iF]->center_mass();
                // _center_mass += _faces[iF]->center_mass();
            }
            _measure /= 3.0;
            _center_mass /= (4.0 * _measure);
            // _center_mass /= 6.0;
        }

        void Cell::set_face_orientations() 
        {
            VectorRd vertex_avg(VectorRd::Zero());
            for(auto & f : _faces)
            {
                vertex_avg += f->center_mass();
            }
            vertex_avg /= _faces.size();
            _face_orientations.reserve(_faces.size());
            for(size_t iF = 0; iF < _faces.size(); ++iF)
            {
                _face_orientations[iF] = Math::sgn(_faces[iF]->normal().dot(_faces[iF]->center_mass() - vertex_avg)); // star-shaped wrt vertex average
            }

        }

        int Cell::face_orientation(const size_t face_index) const
        {
            assert(face_index < _faces.size());
            return (_face_orientations[face_index]);
        }

        VectorRd Cell::face_normal(const size_t face_index) const
        {
            assert(face_index < _faces.size());
            return this->face_orientation(face_index) * (_faces[face_index]->normal());
        }



        bool Cell::test() const
        {
            bool valid = true;
            if(_faces.size() + _vertices.size() != _edges.size() + 2)
            {
                std::cerr << "**** Cell " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") does not satisfy Euler's formula (F + V = E + 2): has " << _vertices.size() << " vertices, " << _edges.size() << " edges, " << _faces.size() << " faces.\n";
                valid = false;
            }
            if(_measure < 1E-14)
            {
                std::cerr << "**** Cell " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has trivial measure: " << _measure << "\n";
                valid = false;
            }

            // integrate nTF over the boundary. Should be zero.
            VectorRd bdry_integral(VectorRd::Zero());
            for(size_t iTF = 0; iTF < _faces.size(); ++iTF)
            {
                bdry_integral += this->face_normal(iTF) * _faces[iTF]->measure();
            }
            if(bdry_integral.norm() > 1E-14)
            {
                std::cerr << "**** Cell " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has ill-formed boundary. Integral of outer normal on boundary has size " << bdry_integral.norm() << "\n";
                valid = false;
            }
            return valid;
        }
    }
}
