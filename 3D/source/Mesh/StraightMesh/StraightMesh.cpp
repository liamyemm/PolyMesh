#include "StraightMesh.hpp"
#include <iostream>

namespace PolyMesh3D
{
    namespace StraightMesh
    {

        Mesh::Mesh() {}
        Mesh::~Mesh()
        {
            for (auto &vertex : _vertices)
            {
                delete vertex;
            }
            for (auto &edge : _edges)
            {
                delete edge;
            }
            for (auto &face : _faces)
            {
                delete face;
            }
            for (auto &cell : _cells)
            {
                delete cell;
            }
        }

        void Mesh::set_name(std::string name) { _mesh_name = name; }
        std::string Mesh::get_name() { return _mesh_name; }

        double Mesh::h_max() const
        {
            double val = 0.0;
            for (auto &cell : _cells)
            {
                val = std::max(val, cell->diam());
            }
            return val;
        }

        std::size_t Mesh::dim() const { return 3; }

        std::size_t Mesh::n_vertices() const { return _vertices.size(); }
        std::size_t Mesh::n_edges() const { return _edges.size(); }
        std::size_t Mesh::n_faces() const { return _faces.size(); }
        std::size_t Mesh::n_cells() const { return _cells.size(); }

        std::size_t Mesh::n_b_vertices() const { return _b_vertices.size(); }
        std::size_t Mesh::n_b_edges() const { return _b_edges.size(); }
        std::size_t Mesh::n_b_faces() const { return _b_faces.size(); }
        std::size_t Mesh::n_b_cells() const { return _b_cells.size(); }

        std::size_t Mesh::n_i_vertices() const { return _i_vertices.size(); }
        std::size_t Mesh::n_i_edges() const { return _i_edges.size(); }
        std::size_t Mesh::n_i_faces() const { return _i_faces.size(); }
        std::size_t Mesh::n_i_cells() const { return _i_cells.size(); }

        std::vector<Vertex *> Mesh::get_vertices() const { return _vertices; }
        std::vector<Edge *> Mesh::get_edges() const { return _edges; }
        std::vector<Face *> Mesh::get_faces() const { return _faces; }
        std::vector<Cell *> Mesh::get_cells() const { return _cells; }

        std::vector<Vertex *> Mesh::get_b_vertices() const { return _b_vertices; }
        std::vector<Edge *> Mesh::get_b_edges() const { return _b_edges; }
        std::vector<Face *> Mesh::get_b_faces() const { return _b_faces; }
        std::vector<Cell *> Mesh::get_b_cells() const { return _b_cells; }

        std::vector<Vertex *> Mesh::get_i_vertices() const { return _i_vertices; }
        std::vector<Edge *> Mesh::get_i_edges() const { return _i_edges; }
        std::vector<Face *> Mesh::get_i_faces() const { return _i_faces; }
        std::vector<Cell *> Mesh::get_i_cells() const { return _i_cells; }

        void Mesh::add_vertex(Vertex *vertex)
        {
            assert(std::find(_vertices.begin(), _vertices.end(), vertex) == _vertices.end());
            _vertices.push_back(vertex);
        }

        void Mesh::add_edge(Edge *edge)
        {
            assert(std::find(_edges.begin(), _edges.end(), edge) == _edges.end());
            _edges.push_back(edge);
        }

        void Mesh::add_face(Face *face)
        {
            assert(std::find(_faces.begin(), _faces.end(), face) == _faces.end());
            _faces.push_back(face);
        }

        void Mesh::add_cell(Cell *cell)
        {
            assert(std::find(_cells.begin(), _cells.end(), cell) == _cells.end());
            _cells.push_back(cell);
        }

        void Mesh::add_b_vertex(Vertex *vertex)
        {
            assert(std::find(_b_vertices.begin(), _b_vertices.end(), vertex) == _b_vertices.end());
            _b_vertices.push_back(vertex);
        }

        void Mesh::add_b_edge(Edge *edge)
        {
            assert(std::find(_b_edges.begin(), _b_edges.end(), edge) == _b_edges.end());
            _b_edges.push_back(edge);
        }

        void Mesh::add_b_face(Face *face)
        {
            assert(std::find(_b_faces.begin(), _b_faces.end(), face) == _b_faces.end());
            _b_faces.push_back(face);
        }

        void Mesh::add_b_cell(Cell *cell)
        {
            assert(std::find(_b_cells.begin(), _b_cells.end(), cell) == _b_cells.end());
            _b_cells.push_back(cell);
        }

        void Mesh::add_i_vertex(Vertex *vertex)
        {
            assert(std::find(_i_vertices.begin(), _i_vertices.end(), vertex) == _i_vertices.end());
            _i_vertices.push_back(vertex);
        }

        void Mesh::add_i_edge(Edge *edge)
        {
            assert(std::find(_i_edges.begin(), _i_edges.end(), edge) == _i_edges.end());
            _i_edges.push_back(edge);
        }

        void Mesh::add_i_face(Face *face)
        {
            assert(std::find(_i_faces.begin(), _i_faces.end(), face) == _i_faces.end());
            _i_faces.push_back(face);
        }

        void Mesh::add_i_cell(Cell *cell)
        {
            assert(std::find(_i_cells.begin(), _i_cells.end(), cell) == _i_cells.end());
            _i_cells.push_back(cell);
        }

        // Note that all these assume that a MeshObject's index is equal to its position in the mesh!!
        Vertex *Mesh::vertex(std::size_t index) const
        {
            assert(index < _vertices.size());
            return _vertices[index];
        }

        Edge *Mesh::edge(std::size_t index) const
        {
            assert(index < _edges.size());
            return _edges[index];
        }

        Face *Mesh::face(std::size_t index) const
        {
            assert(index < _faces.size());
            return _faces[index];
        }

        Cell *Mesh::cell(std::size_t index) const
        {
            assert(index < _cells.size());
            return _cells[index];
        }

        Vertex *Mesh::b_vertex(std::size_t index) const
        {
            assert(index < _b_vertices.size());
            return _b_vertices[index];
        }

        Edge *Mesh::b_edge(std::size_t index) const
        {
            assert(index < _b_edges.size());
            return _b_edges[index];
        }

        Face *Mesh::b_face(std::size_t index) const
        {
            assert(index < _b_faces.size());
            return _b_faces[index];
        }

        Cell *Mesh::b_cell(std::size_t index) const
        {
            assert(index < _b_cells.size());
            return _b_cells[index];
        }

        Vertex *Mesh::i_vertex(std::size_t index) const
        {
            assert(index < _i_vertices.size());
            return _i_vertices[index];
        }

        Edge *Mesh::i_edge(std::size_t index) const
        {
            assert(index < _i_edges.size());
            return _i_edges[index];
        }

        Face *Mesh::i_face(std::size_t index) const
        {
            assert(index < _i_faces.size());
            return _i_faces[index];
        }

        Cell *Mesh::i_cell(std::size_t index) const
        {
            assert(index < _i_cells.size());
            return _i_cells[index];
        }

        void Mesh::renum(const char B, const std::vector<size_t> new_to_old)
        {

            switch (B)
            {
                case 'C':
                {
                    std::vector<Cell *> old_index = _cells;
                    for (size_t i = 0; i < _cells.size(); i++)
                    {
                        old_index[new_to_old[i]]->set_global_index(i);
                        _cells[i] = old_index[new_to_old[i]];
                    }
                    break;
                }

                case 'F':
                {
                    std::vector<Face *> old_index = _faces;
                    for (size_t i = 0; i < _faces.size(); i++)
                    {
                        old_index[new_to_old[i]]->set_global_index(i);
                        _faces[i] = old_index[new_to_old[i]];
                    }
                    break;
                }

                case 'E':
                {
                    std::vector<Edge *> old_index = _edges;
                    for (size_t i = 0; i < _edges.size(); i++)
                    {
                        old_index[new_to_old[i]]->set_global_index(i);
                        _edges[i] = old_index[new_to_old[i]];
                    }
                    break;
                }

                case 'V':
                {
                    std::vector<Vertex *> old_index = _vertices;
                    for (size_t i = 0; i < _vertices.size(); i++)
                    {
                        old_index[new_to_old[i]]->set_global_index(i);
                        _vertices[i] = old_index[new_to_old[i]];
                    }
                    break;
                }
            }
        }

        void Mesh::remove_vertex(Vertex *vertex)
        {
            assert(std::find(_vertices.begin(), _vertices.end(), vertex) != _vertices.end());

            auto pos = std::find(_vertices.begin(), _vertices.end(), vertex);
            _vertices.erase(pos);

            auto pos_i = std::find(_i_vertices.begin(), _i_vertices.end(), vertex);
            if (pos_i != _i_vertices.end())
            {
                _i_vertices.erase(pos_i);
            }

            auto pos_b = std::find(_b_vertices.begin(), _b_vertices.end(), vertex);
            if (pos_b != _b_vertices.end())
            {
                _b_vertices.erase(pos_b);
            }

            for (auto &cell : _cells)
            {
                cell->remove_vertex(vertex);
            }

            for (auto &face : _faces)
            {
                face->remove_vertex(vertex);
            }

            for (auto &edge : _edges)
            {
                edge->remove_vertex(vertex);
            }

            for (auto &vertex : _vertices)
            {
                vertex->remove_vertex(vertex);
            }

            delete vertex;
        }

        void Mesh::remove_edge(Edge *edge)
        {
            assert(std::find(_edges.begin(), _edges.end(), edge) != _edges.end());

            auto pos = std::find(_edges.begin(), _edges.end(), edge);
            _edges.erase(pos);

            auto pos_i = std::find(_i_edges.begin(), _i_edges.end(), edge);
            if (pos_i != _i_edges.end())
            {
                _i_edges.erase(pos_i);
            }

            auto pos_b = std::find(_b_edges.begin(), _b_edges.end(), edge);
            if (pos_b != _b_edges.end())
            {
                _b_edges.erase(pos_b);
            }

            for (auto &cell : _cells)
            {
                cell->remove_edge(edge);
            }

            for (auto &face : _faces)
            {
                face->remove_edge(edge);
            }

            for (auto &vertex : _vertices)
            {
                vertex->remove_edge(edge);
            }

            delete edge;
        }

        void Mesh::remove_face(Face *face)
        {
            assert(std::find(_faces.begin(), _faces.end(), face) != _faces.end());

            auto pos = std::find(_faces.begin(), _faces.end(), face);
            _faces.erase(pos);

            auto pos_i = std::find(_i_faces.begin(), _i_faces.end(), face);
            if (pos_i != _i_faces.end())
            {
                _i_faces.erase(pos_i);
            }

            auto pos_b = std::find(_b_faces.begin(), _b_faces.end(), face);
            if (pos_b != _b_faces.end())
            {
                _b_faces.erase(pos_b);
            }

            for (auto &cell : _cells)
            {
                cell->remove_face(face);
            }

            for (auto &edge : _edges)
            {
                edge->remove_face(face);
            }

            for (auto &vertex : _vertices)
            {
                vertex->remove_face(face);
            }

            delete face;
        }

        void Mesh::remove_cell(Cell *cell)
        {
            assert(std::find(_cells.begin(), _cells.end(), cell) != _cells.end());

            auto pos = std::find(_cells.begin(), _cells.end(), cell);
            _cells.erase(pos);

            auto pos_i = std::find(_i_cells.begin(), _i_cells.end(), cell);
            if (pos_i != _i_cells.end())
            {
                _i_cells.erase(pos_i);
            }

            auto pos_b = std::find(_b_cells.begin(), _b_cells.end(), cell);
            if (pos_b != _b_cells.end())
            {
                _b_cells.erase(pos_b);
            }

            for (auto &face : _faces)
            {
                face->remove_cell(cell);
            }

            for (auto &edge : _edges)
            {
                edge->remove_cell(cell);
            }

            for (auto &vertex : _vertices)
            {
                vertex->remove_cell(cell);
            }

            delete cell;
        }

        double orientation(const VectorRd v1, const VectorRd v2, const VectorRd v3, const VectorRd v4) // returns orientation of tetrahedron v1,v2,v3,v4
        {
            return Math::sgn( (v4 - v1).dot((v2 - v1).cross(v3 - v1)));
        }

        size_t Mesh::find_cell(const VectorRd x) const
        {
            // this algorithm assumes cells and faces are star shaped wrt their respective center masses.

            // Locate neighbouring cells
            std::vector<size_t> neighbouring_cells;
            for (Cell *T : get_cells())
            {
                if ((T->center_mass() - x).norm() <= T->diam())
                {
                    neighbouring_cells.push_back(T->global_index());
                }
            }

            // In neighbouring cells, find one such that x is contained in one of the subtriangles
            size_t i = 0;
            bool found = false;
            size_t iT = 0;
            while (i < neighbouring_cells.size() && !found)
            {
                iT = neighbouring_cells[i];
                Cell *T = cell(iT);
                // the following is a simplicial sub-division of T assuming star-shapedness
                for(auto & F : T->get_faces())
                {
                    for (size_t iV = 0; iV < F->n_vertices(); ++iV)
                    {
                        int o1 = orientation(x, F->vertex(iV)->coords(), F->vertex( (iV + 1) % F->n_vertices() )->coords(), F->center_mass());
                        int o2 = orientation(x, F->vertex( (iV + 1) % F->n_vertices() )->coords(), F->center_mass(), T->center_mass());
                        int o3 = orientation(x, F->center_mass(), T->center_mass(), F->vertex(iV)->coords());
                        int o4 = orientation(x, T->center_mass(), F->vertex(iV)->coords(), F->vertex( (iV + 1) % F->n_vertices() )->coords());

                        found = (o1 >= 0.0 && o2 >= 0.0 && o3 >= 0.0 && o4 >= 0.0) || (o1 <= 0.0 && o2 <= 0.0 && o3 <= 0.0 && o4 <= 0.0); // if all non neg or non pos must be in or on tetrahedron
                        if (found)
                            break;
                    }
                }
                ++i;
            }
            assert(found);

            return iT;
        }


        double Mesh::regularity() const
        {
            // Evaluation of max of ratio "diam of cell / radius ball inscribed in cell"

            double reg = 0.0;
            for (auto &T : _cells)
            {
                double hT = T->diam();
                VectorRd xT = T->center_mass();

                double rhoT = hT;
                for (auto &F : T->get_faces())
                {
                    rhoT = std::min(rhoT, std::abs((xT - F->center_mass()).dot(F->normal()))); // Assuming star-shaped wrt xT
                }
                reg = std::max(reg, hT / rhoT);
            }

            return reg;
        }

        bool Mesh::test() const
        {
            bool vert_valid = true;
            bool edge_valid = true;
            bool face_valid = true;
            bool cell_valid = true;

            std::cout << "\nChecking validity of the mesh\n";

            std::cout << "\n  Checking vertices\n";
            for(size_t iV = 0; iV < _vertices.size(); ++iV)
            {
                vert_valid = vert_valid & _vertices[iV]->test();
                if(_vertices[iV]->global_index() != iV)
                {
                    std::cerr << "    Vertex " << _vertices[iV]->global_index() << " incorrectly indexed.\n";
                    vert_valid = false;
                }
            }
            if(vert_valid)
            {
                std::cout << "    Vertices are OK\n";
            }

            std::cout << "\n  Checking edges\n";
            for(size_t iE = 0; iE < _edges.size(); ++iE)
            {
                edge_valid = edge_valid & _edges[iE]->test();
                if(_edges[iE]->global_index() != iE)
                {
                    std::cerr << "    Edge " << _edges[iE]->global_index() << " incorrectly indexed.\n";
                    edge_valid = false;
                }
            }
            if(edge_valid)
            {
                std::cout << "    Edges are OK\n";
            }

            std::cout << "\n  Checking faces\n";
            for(size_t iF = 0; iF < _faces.size(); ++iF)
            {
                face_valid = face_valid & _faces[iF]->test();
                if(_faces[iF]->global_index() != iF)
                {
                    std::cerr << "    Face " << _faces[iF]->global_index() << " incorrectly indexed.\n";
                    face_valid = false;
                }
            }
            if(face_valid)
            {
                std::cout << "    Faces are OK\n";
            }

            std::cout << "\n  Checking cells\n";
            for(size_t iT = 0; iT < _cells.size(); ++iT)
            {
                cell_valid = cell_valid & _cells[iT]->test();
                if(_cells[iT]->global_index() != iT)
                {
                    std::cerr << "    Cell " << _cells[iT]->global_index() << " incorrectly indexed.\n";
                    cell_valid = false;
                }
            }
            if(cell_valid)
            {
                std::cout << "    Cells are OK\n";
            }

            bool valid = (vert_valid && edge_valid && face_valid && cell_valid);
            if(vert_valid && edge_valid && face_valid && cell_valid)
            {
                std::cout << "\n  Mesh is OK\n";
            }
            else
            {
                std::cout << "\n  Mesh is not OK\n";
            }
            return valid;
        }
    }
}