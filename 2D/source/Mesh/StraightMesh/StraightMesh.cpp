#include "StraightMesh.hpp"
#include <iostream>

namespace PolyMesh2D
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

        void Mesh::plot_mesh(std::ofstream *out) const
        {
            for (auto &edge : _edges)
            {
                edge->plot(out);
            }
        }

        void Mesh::plot_triangulation(std::ofstream *out) const
        {
            for (auto &cell : _cells)
            {
                cell->plot_triangulation(out);
            }
        }

        void Mesh::remove_edge(Edge *edge)
        {
            assert(edge->n_cells() == 0); // can only remove edges with no cells
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

            for (auto &vertex : _vertices)
            {
                vertex->remove_edge(edge);
            }

            delete edge;
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

        size_t Mesh::find_cell(const VectorRd x) const
        {

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

                for (auto &simplex : T->triangulation())
                {
                    double area1 = signed_area(x, simplex[0], simplex[1]);
                    double area2 = signed_area(x, simplex[1], simplex[2]);
                    double area3 = signed_area(x, simplex[2], simplex[0]);

                    found = (area1 >= 0.0 && area2 >= 0.0 && area3 >= 0.0) || (area1 <= 0.0 && area2 <= 0.0 && area3 <= 0.0); // if all non neg or non pos must be in or on triangle
                    if (found)
                        break;
                }
                ++i;
            }
            assert(found);

            return iT;
        }


        std::array<double, 2> Mesh::regularity()
        {
            /// Regularity factor =
            ///   1st component: maximum of
            ///      * diameter of cell / (measure of cell)^{1/dim}
            ///      * diameter of cell / diameter of edge  [for each edge of the cell]
            ///
            ///   2nd component: evaluation of max of ratio "diam of cell / radius ball inscribed in cell"

            std::vector<std::vector<double>> reg_cell(n_cells(), {0.0, 0.0});
            std::size_t count = 0;
            for (auto &T : _cells)
            {
                double hT = T->diam();
                VectorRd xT = T->center_mass();

                reg_cell[count][0] = hT / pow(T->measure(), 1.0 / 2.0);

                double rhoT = hT;
                std::vector<Edge *> edges = T->get_edges();
                for (auto &F : edges)
                {
                    double hF = F->diam();
                    VectorRd xF = F->center_mass();
                    VectorRd nTF = F->normal(); // sign does not matter

                    reg_cell[count][0] = std::max(reg_cell[count][0], hT / hF);

                    rhoT = std::min(rhoT, std::abs((xT - xF).dot(nTF))); // If xT is not in T, is this really a good measure?
                }
                reg_cell[count][1] = hT / rhoT;
                ++count; // could just use iterators
            }

            std::array<double, 2> value{2, 0.0};
            for (size_t iT = 0; iT < n_cells(); iT++)
            {
                value[0] = std::max(value[0], reg_cell[iT][0]);
                value[1] = std::max(value[1], reg_cell[iT][1]);
            }

            return value;
        }
    }
}