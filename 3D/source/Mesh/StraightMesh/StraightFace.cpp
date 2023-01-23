#include <iostream>

#include "StraightVertex.hpp"
#include "StraightEdge.hpp"
#include "StraightFace.hpp"
#include "StraightCell.hpp"

namespace PolyMesh3D
{
    namespace StraightMesh
    {
        Face::Face(size_t index, std::vector<Edge *> edges)
            : Polytope::Polytope(index)
        {
            _edges = edges;
            // Compute the normal using the most perpendicular pair of edges.
            double norm = 0.0;
            VectorRd perp_vec = VectorRd::Zero();
            for (auto &e1 : _edges)
            {
                for (auto &e2 : _edges)
                {
                    if (e1 == e2)
                    {
                        continue;
                    }
                    VectorRd vec1(e1->tangent());
                    VectorRd vec2(e2->tangent());
                    if ((vec1.cross(vec2)).norm() > norm)
                    {
                        perp_vec = vec1.cross(vec2);
                        norm = perp_vec.norm();
                    }
                }
            }
            _normal = perp_vec / norm;

            // Set outer normals. Will initially either be all outer, or all inner. If all inner, _measure will be negative, so switch them.
            // Following algorithm assumes _edges is ordered (be it clock-wise or counter-clockwise wrt the normal). If _edges s randomly ordered, need to add an ordering routine

            _outer_normals.reserve(_edges.size()); // VectorRd is a static type, so memory requirement known at compile time.
            bool prior_was_switched = false; 
            for (size_t iE = 0; iE < _edges.size(); ++iE)
            {
                VectorRd tan_vec = (_edges[iE]->coords()[1] - _edges[iE]->coords()[0]).normalized();
                if(iE > 0)
                {
                    // test if tan_vec is oriented the same as previous tan_vec. Otherwise switch it.
                    if( (_edges[iE]->coords()[0] != _edges[iE - 1]->coords()[1]) && (_edges[iE]->coords()[1] != _edges[iE - 1]->coords()[0]) )
                    {
                        assert( (_edges[iE]->coords()[0] == _edges[iE - 1]->coords()[0]) || (_edges[iE]->coords()[1] == _edges[iE - 1]->coords()[1]) );
                        if(!prior_was_switched)
                        {
                            tan_vec *= -1.0;
                            prior_was_switched = true;
                        }
                        else
                        {
                            prior_was_switched = false;
                        }
                    }
                    else
                    {
                        if(prior_was_switched)
                        {
                            tan_vec *= -1.0;
                            prior_was_switched = true;
                        }
                        else
                        {
                            prior_was_switched = false;
                        }
                    }
                }
                _outer_normals.push_back(tan_vec.cross(_normal));
                _measure += _edges[iE]->measure() * _outer_normals[iE].dot(_edges[iE]->center_mass()); // constant integral on the edge so can be evaluated at edge center

                // while we are looping over edges, may as well find diameter
                for (size_t jE = iE + 1; jE < _edges.size(); ++jE)
                {
                    _diameter = std::max(_diameter, (_edges[jE]->coords()[0] - _edges[iE]->coords()[0]).norm());
                    // vertices in each edge may not be listed the same way, so need the following to guarentee capturing every vertex pair
                    // Not the most efficient, but oh well...
                    _diameter = std::max(_diameter, (_edges[jE]->coords()[0] - _edges[iE]->coords()[1]).norm());
                }
            }
            _measure /= 2.0; // Faces are 2D, so dividing by dimension.

            double multiplier = (_measure < 0.0 ? -1.0 : 1.0);
            _measure *= multiplier;
            // _outer_normals[0] *= multiplier;

            // if _measure is zero (or sufficiently close to it) then there will be many issues. Tested for in test()

            // find the center mass and invert outer_normals if facing inwards
            _center_mass = VectorRd::Zero();
            for (size_t iE = 0; iE < _edges.size(); ++iE)
            {
                _outer_normals[iE] *= multiplier;

                double x_minus_x0_dot_n = (_edges[iE]->center_mass() - _edges[0]->center_mass()).dot(_outer_normals[iE]);
                VectorRd point = (1.0 / 3.0) * _edges[iE]->center_mass() + (1.0 / 6.0) * _edges[0]->center_mass();

                _center_mass += _edges[iE]->measure() * x_minus_x0_dot_n * point;
            }
            _center_mass /= _measure;

            set_edge_orientations(); // set the edge orientations
        }

        void Face::set_edge_orientations()
        {
            _edge_orientations.clear();                // reset to empty vector
            _edge_orientations.reserve(_edges.size()); // VectorRd is a static type, so memory requirement known at compile time.

            for (size_t iE = 0; iE < _edges.size(); ++iE)
            {
                VectorRd normal = this->edge_normal(iE);
                _edge_orientations.push_back(Math::sgn(normal.dot(_outer_normals[iE])));
            }

            // _outer_normals.clear(); // no longer needed
        }

        int Face::edge_orientation(const size_t edge_index) const
        {
            assert(edge_index < _edges.size());
            return (_edge_orientations[edge_index]);
        }

        VectorRd Face::edge_normal(const size_t edge_index) const
        {
            assert(edge_index < _edges.size());
            return _normal.cross(_edges[edge_index]->tangent());
        }

        VectorRd Face::normal() const
        {
            return _normal;
        }

        bool Face::test() const
        {
            bool valid = true;
            if(_vertices.size() != _edges.size())
            {
                std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has " << _vertices.size() << " vertices and" << _edges.size() << " edges.\n";
                valid = false;
            }
            if(_measure < 1e-12)
            {
                std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has trivial measure: " << _measure << "\n";
                valid = false;
            }

            // integrate nFE over the boundary of face. Should be zero.
            VectorRd bdry_integral(VectorRd::Zero());
            for(size_t iFE = 0; iFE < _edges.size(); ++iFE)
            {
                bdry_integral += this->edge_orientation(iFE) * this->edge_normal(iFE) * _edges[iFE]->measure();
            }
            if(bdry_integral.norm() > 1e-12)
            {
                std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has ill-formed boundary. Integral of outer normal on boundary has size " << bdry_integral.norm() << "\n";
                valid = false;
            }
            if (is_boundary())
            {
                if (_cells.size() != 1)
                {
                    std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") is a boundary face but has " << _cells.size() << " cells.\n";
                    valid = false;
                }
            }
            else
            {
                if (_cells.size() != 2)
                {
                    std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") is an internal face but has " << _cells.size() << " cells.\n";
                    valid = false;
                }
                else
                {
                    // Get orientations of adjoining cells and check they multiply to one.
                    int index_1 = -1;
                    int index_2 = -1;

                    for(size_t iTF = 0; iTF < this->cell(0)->n_faces(); ++iTF)
                    {
                        if(this->cell(0)->face(iTF) == this)
                        {
                            index_1 = iTF;
                            break;
                        }
                    }
                    for(size_t iTF = 0; iTF < this->cell(1)->n_faces(); ++iTF)
                    {
                        if(this->cell(1)->face(iTF) == this)
                        {
                            index_2 = iTF;
                            break;
                        }
                    }

                    assert(index_1 >= 0 && index_2 >= 0);

                    // test if internal normals are oriented correctly

                    if(this->cell(0)->face_orientation(index_1) * this->cell(1)->face_orientation(index_2) != -1)
                    {
                        std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") has orientation issues.\n";
                        std::cerr << "   cell/center/orientation: " << std::endl;
                        std::cerr << "      " << this->cell(0)->global_index() << " / (" << this->cell(0)->center_mass().transpose() << ") / " << this->cell(0)->face_orientation(index_1) << std::endl;
                        std::cerr << "      " << this->cell(1)->global_index() << " / (" << this->cell(1)->center_mass().transpose() << ") / " << this->cell(1)->face_orientation(index_2) << std::endl << std::endl;
                        valid = false;
                    }
                }
            }
            
            // Checking that the face is planar: all vectors in the face should be perpendicular to the normal
            for (size_t iE = 0; iE < _edges.size(); iE++)
            {
                VectorRd tan_vec = _edges[iE]->tangent();
                VectorRd center_vec = (_edges[iE]->center_mass() - _center_mass).normalized();
                if (std::abs(tan_vec.dot(_normal)) > 1e-12 || std::abs(center_vec.dot(_normal)) > 1e-8)
                {
                    std::cerr << "**** Face " << _index << " with center at (" << _center_mass(0) << ", " << _center_mass(1) << ", " << _center_mass(2) << ") is non-planar.\n";

                    std::cerr << "******** tan.normal: " << std::abs(tan_vec.dot(_normal)) << " center.normal: " << std::abs(center_vec.dot(_normal)) << "\n";
                    valid = false;
                    break;
                }
            }

            // if(!valid)
            // {
            //         std::cout << "\n";
            //         for(auto & v : _vertices)
            //         {
            //             std::cout << v->coords().transpose() << "\n";
            //         }
            //         std::cout << "\n";
            //         // std::cout << _normal.transpose() << "\n";
            //         // std::cout << "\n";
            //         for(size_t iE = 0; iE < _edges.size(); ++iE)
            //         {
            //             std::cout << _edges[iE]->center_mass().transpose() << "\n";
            //             std::cout << _outer_normals[iE].transpose() << "\n";
            //             std::cout << "\n";
            //         }
            //         std::cout << "\n";
            // }
            return valid;
        }
    }
}
