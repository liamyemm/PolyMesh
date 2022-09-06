#include <../Mesh2D/Polytope2D.hpp>

#ifndef _FUNCTION2D_HPP
#define _FUNCTION2D_HPP

namespace Mesh2D
{
    template <size_t object_dim>
    class CurvedObject
    {
    public:
        ///< Constructor for a Polytope defined by its simplices
        CurvedObject(size_t index);

        ///< Null constructor
        CurvedObject();

        ///< Destructor
        ~CurvedObject();

        // inline size_t dim() const { return object_dim; }

        inline size_t global_index() const { return _index; }                     ///< Return the global index of the Polytope
        inline double diam() const { return _diameter; }                          ///< Return the diameter of the Polytope
        inline VectorRd center_mass() const { return _center_mass; }              ///< Return the center mass of the Polytope
        inline double measure() const { return _measure; }                        ///< Return the Lebesgue measure of the Polytope
        inline Simplices<object_dim> get_simplices() const { return _simplices; } ///< Return the simplices making up the Polytope
        inline void set_global_index(const size_t idx) { _index = idx; }          ///< Set the global index
        inline bool is_boundary() const { return _is_boundary; }                  ///< Return true if Polytope is a boundary object, false otherwise
        inline void set_boundary(bool val) { _is_boundary = val; }                ///< Set the boundary value of the Polytope

        inline std::vector<Polytope<0> *> get_vertices() const { return _vertices; }       ///< Return the vertices of the Polytope
        inline std::vector<Polytope<1> *> get_edges() const { return _edges; }             ///< Return the edges of the Polytope
        inline std::vector<Polytope<DIMENSION - 1> *> get_faces() const { return _edges; } ///< Return the faces of the Polytope
        inline std::vector<Polytope<DIMENSION> *> get_cells() const { return _cells; }     ///< Return the cells of the Polytope

        inline size_t n_vertices() const { return _vertices.size(); } ///< Return the number of vertices of the Polytope
        inline size_t n_edges() const { return _edges.size(); }       ///< Return the number of edges of the Polytope
        inline size_t n_faces() const { return _edges.size(); }       ///< Return the number of faces of the Polytope
        inline size_t n_cells() const { return _cells.size(); }       ///< Return the number of cells of the Polytope

        Polytope<0> *vertex(const size_t i) const;           ///< Return the i-th vertex of the Polytope
        Polytope<1> *edge(const size_t i) const;             ///< Return the i-th edge of the Polytope
        Polytope<DIMENSION - 1> *face(const size_t i) const; ///< Return the i-th face of the Polytope
        Polytope<DIMENSION> *cell(const size_t i) const;     ///< Return the i-th cell of the Polytope

        void add_vertex(Polytope<0> *vertex);         ///< Add a vertex to the Polytope
        void add_edge(Polytope<1> *edge);             ///< Add an edge to the Polytope
        void add_face(Polytope<DIMENSION - 1> *face); ///< Add a face to the Polytope
        void add_cell(Polytope<DIMENSION> *cell);     ///< Add a cell to the Polytope

        int index_vertex(const Polytope<0> *vertex) const;         ///< Returns the local index of a vertex
        int index_edge(const Polytope<1> *edge) const;             ///< Returns the local index of an edge
        int index_face(const Polytope<DIMENSION - 1> *face) const; ///< Returns the local index of a face
        int index_cell(const Polytope<DIMENSION> *cell) const;     ///< Returns the local index of a cell

        VectorRd coords() const; ///< Return the coordinates of a Vertex

        VectorRd face_normal(const size_t face_index) const; ///< Return the outer normal of a Cell towards the Face located at face_index
        VectorRd edge_normal(const size_t edge_index) const; ///< Return the edge normal of a 2D object

        int face_orientation(const size_t face_index) const; ///< Return the orientation of a Face
        int edge_orientation(const size_t edge_index) const; ///< Return the orientation of a Edge
        int vertex_orientation(const size_t vertex_index) const; ///< Return the orientation of a Vertex

        VectorRd normal() const;  ///< Return the normal of a Face
        VectorRd tangent() const; ///< Return the tangent of a Edge

        void construct_face_normals(); ///< Set the directions of the face normals of a cell

        void plot_simplices(std::ofstream *out) const; ///< Plot the simplices to out
        void plot(std::ofstream *out) const; ///< Plot the polytope to out

    private:
        size_t _index;
        VectorRd _center_mass;
        double _measure;
        double _diameter;
        bool _is_boundary;
        Simplices<object_dim> _simplices;
        VectorRd _normal;                  // uninitialised unless object_dim == DIMENSION - 1 (face)
        VectorRd _tangent;                  // uninitialised unless object_dim == 1 (edge)
        std::vector<int> _face_directions; // empty unless object_dim == DIMENSION (cell)

        std::vector<Polytope<0> *> _vertices;
        std::vector<Polytope<1> *> _edges;
        // std::vector<Polytope<DIMENSION - 1>*> _faces;
        std::vector<Polytope<DIMENSION> *> _cells;
    };
}

#endif