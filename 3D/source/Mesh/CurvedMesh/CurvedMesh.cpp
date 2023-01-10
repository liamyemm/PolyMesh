#include "CurvedMesh.hpp"
#include <iostream>
#include "GaussLegendre.hpp"

namespace PolyMesh2D
{
    namespace CurvedMesh
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

        std::size_t Mesh::dim() const { return 2; }

        std::size_t Mesh::n_vertices() const { return _vertices.size(); }
        std::size_t Mesh::n_edges() const { return _edges.size(); }
        std::size_t Mesh::n_cells() const { return _cells.size(); }

        std::size_t Mesh::n_b_vertices() const { return _b_vertices.size(); }
        std::size_t Mesh::n_b_edges() const { return _b_edges.size(); }
        std::size_t Mesh::n_b_cells() const { return _b_cells.size(); }

        std::size_t Mesh::n_i_vertices() const { return _i_vertices.size(); }
        std::size_t Mesh::n_i_edges() const { return _i_edges.size(); }
        std::size_t Mesh::n_i_cells() const { return _i_cells.size(); }

        std::vector<Vertex *> Mesh::get_vertices() const { return _vertices; }
        std::vector<Edge *> Mesh::get_edges() const { return _edges; }
        std::vector<Cell *> Mesh::get_cells() const { return _cells; }

        std::vector<Vertex *> Mesh::get_b_vertices() const { return _b_vertices; }
        std::vector<Edge *> Mesh::get_b_edges() const { return _b_edges; }
        std::vector<Cell *> Mesh::get_b_cells() const { return _b_cells; }

        std::vector<Vertex *> Mesh::get_i_vertices() const { return _i_vertices; }
        std::vector<Edge *> Mesh::get_i_edges() const { return _i_edges; }
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

            case 'E':
            {
                std::vector<Edge *> old_index = _edges;
                for (size_t i = 0; i < _edges.size(); i++)
                {
                    // if(_edges[i] == nullptr)
                    // {
                    //     std::cout << "gotcha\n";
                    //     exit(1);
                    // }
                    if(old_index[new_to_old[i]] == nullptr)
                    {
                        std::cout << "gotcha\n";
                        exit(1);
                    }
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

        void Mesh::plot_mesh(std::ofstream *out, int partitions) const
        {
            for (auto &edge : _edges)
            {
                edge->plot(out, partitions);
            }
        }

        void Mesh::remove_edge(Edge *edge)
        {
            // assert(edge->n_cells() == 0); // can only remove edges with no cells
            assert(std::find(_edges.begin(), _edges.end(), edge) != _edges.end());

            size_t current_index = edge->global_index();
            for(size_t id = current_index + 1; id < _edges.size(); ++id)
            {
                _edges[id]->set_global_index(id - 1);
            }

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

            for (auto &vertex : edge->get_vertices())
            {
                vertex->remove_edge(edge);
            }

            for (auto &cell : edge->get_cells())
            {
                cell->remove_edge(edge);
            }

            edge->vertex(0)->remove_vertex(edge->vertex(1));
            edge->vertex(1)->remove_vertex(edge->vertex(0));

            // for (auto &cell : edge->get_cells())
            // {
            //     cell->remove_edge(edge);
            // }

            delete edge;
        }

        void Mesh::remove_cell(Cell *cell)
        {
            assert(std::find(_cells.begin(), _cells.end(), cell) != _cells.end());

            size_t current_index = cell->global_index();
            for(size_t id = current_index + 1; id < _cells.size(); ++id)
            {
                _cells[id]->set_global_index(id - 1);
            }

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

            for (auto &edge : cell->get_edges())
            {
                edge->remove_cell(cell);
            }

            for (auto &vertex : cell->get_vertices())
            {
                vertex->remove_cell(cell);
            }

            delete cell;
        }


        bool Mesh::test() const
        {
            bool valid = true;

            if(_i_vertices.size() + _b_vertices.size() != _vertices.size())
            {
                std::cerr << "Error! There are " << _i_vertices.size() << " internal vertices, " << _b_vertices.size() << " boundary vertices and " << _vertices.size() << " total vertices.\n";
                valid = false;
            }

            if(_i_edges.size() + _b_edges.size() != _edges.size())
            {
                std::cerr << "Error! There are " << _i_edges.size() << " internal edges, " << _b_edges.size() << " boundary edges and " << _edges.size() << " total edges.\n";
                valid = false;
            }

            if(_i_cells.size() + _b_cells.size() != _cells.size())
            {
                std::cerr << "Error! There are " << _i_cells.size() << " internal cells, " << _b_cells.size() << " boundary cells and " << _cells.size() << " total cells.\n";
                valid = false;
            }

            for(size_t iE = 0; iE < _i_edges.size(); ++iE)
            {
                int index_1 = -1;
                int index_2 = -1;

                for(size_t iTE = 0; iTE < _i_edges[iE]->cell(0)->n_edges(); ++iTE)
                {
                    if(_i_edges[iE]->cell(0)->edge(iTE) == _i_edges[iE])
                    {
                        index_1 = iTE;
                        break;
                    }
                }
                for(size_t iTE = 0; iTE < _i_edges[iE]->cell(1)->n_edges(); ++iTE)
                {
                    if(_i_edges[iE]->cell(1)->edge(iTE) == _i_edges[iE])
                    {
                        index_2 = iTE;
                        break;
                    }
                }

                assert(index_1 >= 0 && index_2 >= 0);

                // test if internal normals are oriented correctly

                if(_i_edges[iE]->cell(0)->edge_orientation(index_1) * _i_edges[iE]->cell(1)->edge_orientation(index_2) != -1)
                {
                    std::cerr << "Error! Edge " << _i_edges[iE]->global_index() << " has normal issues.\n";
                    valid = false;
                }
            }

            for(size_t iV = 0; iV < _vertices.size(); ++iV)
            {
                valid = valid & _vertices[iV]->test();
                if(_vertices[iV]->global_index() != iV)
                {
                    std::cerr << "Error! Vertex " << _vertices[iV]->global_index() << " incorrectly indexed.\n";
                    valid = false;
                }
            }
            for(size_t iE = 0; iE < _edges.size(); ++iE)
            {
                valid = valid & _edges[iE]->test();
                if(_edges[iE]->global_index() != iE)
                {
                    std::cerr << "Error! Edge " << _edges[iE]->global_index() << " incorrectly indexed.\n";
                    valid = false;
                }
            }
            for(size_t iT = 0; iT < _cells.size(); ++iT)
            {
                valid = valid & _cells[iT]->test();
                if(_cells[iT]->global_index() != iT)
                {
                    std::cerr << "Error! Cell " << _cells[iT]->global_index() << " incorrectly indexed.\n";
                    valid = false;
                }
            }

            // Ensure boundary normals are oriented correctly by testing that \sum_{F\in\Fhb} \int_F v . nTF = 0 for a constant vector v

            double integral = 0.0;
            Eigen::Vector2d const_vec(1.0, 1.0);

            unsigned quad_degree = 20;

            Quadrature::GaussLegendre1D quad(quad_degree);

            for (size_t ibE = 0; ibE < _b_edges.size(); ++ibE)
            {
                Functional::Curve edge_param = _b_edges[ibE]->parameterisation();

                double a = edge_param.tmin;
                double b = edge_param.tmax;

                double edge_len = b - a;
                assert(edge_len > 0.0);

                int normal_direction;

                for(size_t iTE = 0; iTE < _b_edges[ibE]->cell(0)->n_edges(); ++iTE)
                {
                    if(_b_edges[ibE]->cell(0)->edge(iTE) == _b_edges[ibE])
                    {
                        normal_direction = _b_edges[ibE]->cell(0)->edge_orientation(iTE);
                        break;
                    }
                }

                for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
                {
                    double point = edge_len * quad.point(iqn) + a;
                    integral += normal_direction * edge_len * quad.weight(iqn) * const_vec.dot(_b_edges[ibE]->normal(point)) * edge_param.derivative(point).norm();
                }
            }

            if (std::abs(integral) > 1E-12)
            {
                std::cout << "Error! Mesh has boundary integral evaluated to " << integral << ".\n";
                valid = false;
            }
            return valid;
        }
    }
}