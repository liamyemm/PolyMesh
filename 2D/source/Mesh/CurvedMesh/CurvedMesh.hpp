#include "CurvedPolytope.hpp"
#include "CurvedVertex.hpp"
#include "CurvedEdge.hpp"
#include "CurvedCell.hpp"

#ifndef _CURVEDMESH_HPP
#define _CURVEDMESH_HPP

namespace PolyMesh2D
{
    namespace CurvedMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Mesh
        {
        public:
            Mesh();  ///< Null constructor.
            ~Mesh(); ///< Destructor ensuring the deletion of the dynamically allocated vertices, edges, and cells.

            Mesh(const Mesh &) = delete;              ///< Mesh should not be copied
            Mesh(Mesh &&) noexcept = delete;          ///< Mesh should not be moved
            Mesh &operator = (const Mesh &) = delete;     ///< Mesh should not be copied
            Mesh &operator = (Mesh &&) noexcept = delete; ///< Mesh should not be moved

            void set_name(std::string name); ///< set the name of the mesh
            std::string get_name();          ///< getter for the mesh name

            double h_max() const; ///< max diameter of cells

            std::size_t dim() const; ///< dimension of the mesh

            std::size_t n_vertices() const; ///< number of vertices in the mesh.
            std::size_t n_edges() const;    ///< number of edges in the mesh.
            std::size_t n_cells() const;    ///< number of cells in the mesh.

            std::size_t n_b_vertices() const; ///< number of boundary vertices in the mesh.
            std::size_t n_b_edges() const;    ///< number of boundary edges in the mesh.
            std::size_t n_b_cells() const;    ///< number of boundary cells in the mesh.

            std::size_t n_i_vertices() const; ///< number of internal vertices in the mesh.
            std::size_t n_i_edges() const;    ///< number of internal edges in the mesh.
            std::size_t n_i_cells() const;    ///< number of internal cells in the mesh.

            std::vector<Vertex *> get_vertices() const; ///< lists the vertices in the mesh.
            std::vector<Edge *> get_edges() const;      ///< lists the edges in the mesh.
            std::vector<Cell *> get_cells() const;      ///< lists the cells in the mesh.

            std::vector<Vertex *> get_b_vertices() const; ///< lists the boundary vertices in the mesh.
            std::vector<Edge *> get_b_edges() const;      ///< lists the boundary edges in the mesh.
            std::vector<Cell *> get_b_cells() const;      ///< lists the boundary cells in the mesh.

            std::vector<Vertex *> get_i_vertices() const; ///< lists the internal vertices in the mesh.
            std::vector<Edge *> get_i_edges() const;      ///< lists the internal edges in the mesh.
            std::vector<Cell *> get_i_cells() const;      ///< lists the internal cells in the mesh.

            void add_vertex(Vertex *vertex); ///<  adds a vertex to the mesh
            void add_edge(Edge *edge);       ///<  adds a edge to the mesh
            void add_cell(Cell *cell);       ///<  adds a cell to the mesh

            void remove_vertex(Vertex *vertex); ///<  remove a vertex from the mesh
            void remove_edge(Edge *edge); ///<  remove a edge from the mesh
            void remove_cell(Cell *cell); ///<  remove a cell from the mesh

            void add_b_vertex(Vertex *vertex); ///<  adds a boundary vertex to the mesh
            void add_b_edge(Edge *edge);       ///<  adds a boundary edge to the mesh
            void add_b_cell(Cell *cell);       ///<  adds a boundary cell to the mesh

            void add_i_vertex(Vertex *vertex); ///<  adds an internal vertex to the mesh
            void add_i_edge(Edge *edge);       ///<  adds an internal edge to the mesh
            void add_i_cell(Cell *cell);       ///<  adds an internal cell to the mesh

            Vertex *vertex(std::size_t index) const; ///<  get a constant pointer to a vertex using its global index
            Edge *edge(std::size_t index) const;     ///<  get a constant pointer to a edge using its global index
            Cell *cell(std::size_t index) const;     ///<  get a constant pointer to a cell using its global index

            Vertex *b_vertex(std::size_t index) const; ///<  get a constant pointer to a boundary vertex using an index
            Edge *b_edge(std::size_t index) const;     ///<  get a constant pointer to boundary a edge using an index
            Cell *b_cell(std::size_t index) const;     ///<  get a constant pointer to boundary a cell using an index

            Vertex *i_vertex(std::size_t index) const; ///<  get a constant pointer to an internal vertex using an index
            Edge *i_edge(std::size_t index) const;     ///<  get a constant pointer to an internal edge using an index
            Cell *i_cell(std::size_t index) const;     ///<  get a constant pointer to an internal cell using an index

            void renum(const char B, const std::vector<size_t> new_to_old);

            void plot_mesh(std::ofstream *out, int partitions) const;

            bool test() const;

        private:
            std::string _mesh_name;

            std::vector<Vertex *> _vertices;
            std::vector<Edge *> _edges;
            std::vector<Cell *> _cells;

            std::vector<Vertex *> _b_vertices;
            std::vector<Edge *> _b_edges;
            std::vector<Cell *> _b_cells;

            std::vector<Vertex *> _i_vertices;
            std::vector<Edge *> _i_edges;
            std::vector<Cell *> _i_cells;
        };

        //@}
    }
}

#endif
