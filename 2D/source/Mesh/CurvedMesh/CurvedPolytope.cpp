#include "CurvedPolytope.hpp"

namespace PolyMesh2D
{
    namespace CurvedMesh
    {
        Polytope::Polytope(size_t index) : _index(index), _center_mass(VectorRd::Zero()), _measure(0.0), _diameter(0.0), _is_boundary(false) {}
        size_t Polytope::global_index() const { return _index; }
        double Polytope::diam() const { return _diameter; }
        VectorRd Polytope::center_mass() const { return _center_mass; }
        double Polytope::measure() const { return _measure; }
        void Polytope::set_global_index(const size_t idx) { _index = idx; }
        bool Polytope::is_boundary() const { return _is_boundary; }
        void Polytope::set_boundary(bool val) { _is_boundary = val; }

        std::vector<Vertex *> Polytope::get_vertices() const { return _vertices; }
        std::vector<Edge *> Polytope::get_edges() const { return _edges; }
        std::vector<Cell *> Polytope::get_cells() const { return _cells; }

        size_t Polytope::n_vertices() const { return _vertices.size(); }
        size_t Polytope::n_edges() const { return _edges.size(); }
        size_t Polytope::n_cells() const { return _cells.size(); }

        void Polytope::add_vertex(Vertex *vertex)
        {
            assert(std::find(_vertices.begin(), _vertices.end(), vertex) == _vertices.end()); // ensure vertex does not already exist in _vertices
            _vertices.push_back(vertex);
        }

        void Polytope::add_edge(Edge *edge)
        {
            assert(std::find(_edges.begin(), _edges.end(), edge) == _edges.end());
            _edges.push_back(edge);
        }

        void Polytope::add_cell(Cell *cell)
        {
            assert(std::find(_cells.begin(), _cells.end(), cell) == _cells.end());
            _cells.push_back(cell);
        }

        void Polytope::remove_vertex(Vertex *vertex)
        {
            auto pos = std::find(_vertices.begin(), _vertices.end(), vertex);
            if (pos != _vertices.end())
            {
                _vertices.erase(pos);
            }
        }

        void Polytope::remove_edge(Edge *edge)
        {
            auto pos = std::find(_edges.begin(), _edges.end(), edge);
            if (pos != _edges.end())
            {
                _edges.erase(pos);
            }
        }

        void Polytope::remove_cell(Cell *cell)
        {
            auto pos = std::find(_cells.begin(), _cells.end(), cell);
            if (pos != _cells.end())
            {
                _cells.erase(pos);
            }
        }

        Vertex *Polytope::vertex(const size_t i) const
        {
            assert(i < _vertices.size());
            return _vertices[i];
        }

        Edge *Polytope::edge(const size_t i) const
        {
            assert(i < _edges.size());
            return _edges[i];
        }

        Cell *Polytope::cell(const size_t i) const
        {
            assert(i < _cells.size());
            return _cells[i];
        }

        int Polytope::index_vertex(const Vertex *vertex) const
        {
            auto itr = std::find(_vertices.begin(), _vertices.end(), vertex);
            if (itr != _vertices.end())
            {
                return itr - _vertices.begin();
            }
            else
            {
                throw "Vertex not found";
            }
        }

        int Polytope::index_edge(const Edge *edge) const
        {
            auto itr = std::find(_edges.begin(), _edges.end(), edge);
            if (itr != _edges.end())
            {
                return itr - _edges.begin();
            }
            else
            {
                throw "Edge not found";
            }
        }

        int Polytope::index_cell(const Cell *cell) const
        {
            auto itr = std::find(_cells.begin(), _cells.end(), cell);
            if (itr != _cells.end())
            {
                return itr - _cells.begin();
            }
            else
            {
                throw "Cell not found";
            }
        }
    }
}