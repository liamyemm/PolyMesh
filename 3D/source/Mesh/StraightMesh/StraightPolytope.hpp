#include "math.hpp"

// standard libraries
#include <vector>    // std::vector
#include <algorithm> // std::find
#include <fstream>   // std::ofstream

// external linking to Eigen is required
#include <Eigen/Dense> // Eigen::Matrix, Eigen::MatrixXd, Eigen::FullPivLU

#ifndef _POLYTOPE_HPP
#define _POLYTOPE_HPP

/*!
 * @defgroup Mesh
 * @brief Classes and methods for all objects (Vertex, Edge, Face, Cell, Mesh) in a Mesh.
 */

namespace PolyMesh3D
{
    /** The namespace for all objects in a straight mesh. **/
    namespace StraightMesh
    {
        class Vertex;
        class Edge;
        class Face;
        class Cell;

        /*!
         * \addtogroup Mesh
         * @{
         */

        using VectorRd = Eigen::Vector3d; ///< Alias for a 3d vector. Used for coordinates.

        // // ----------------------------------------------------------------------------
        // //                            Class definition
        // // ----------------------------------------------------------------------------

        /**
         * Polytope is the base class for each of Vertex, Edge, Face and Cell and contains methods that are common between all mesh objects
         * Polytope should never be used directly, only via its inherited classes Vertex, Edge, Face and Cell.
         **/
        class Polytope
        {
        public:
            ///< Constructor for the base class. Sets _index to @arg index, _center_mass, _measure, _diameter to zero, and _is_boundary to false.
            Polytope(size_t index);
            virtual ~Polytope();

            size_t global_index() const;             ///< Return the global index of the StraightObject
            double diam() const;                     ///< Return the diameter of the StraightObject
            VectorRd center_mass() const;            ///< Return the center mass of the StraightObject
            double measure() const;                  ///< Return the Lebesgue measure of the StraightObject
            void set_global_index(const size_t idx); ///< Set the global index
            bool is_boundary() const;                ///< Return true if StraightObject is a boundary object, false otherwise
            void set_boundary(bool val);             ///< Set the boundary value of the StraightObject

            std::vector<Vertex *> get_vertices() const; ///< Return the vertices of the StraightObject
            std::vector<Edge *> get_edges() const;      ///< Return the edges of the StraightObject
            std::vector<Face *> get_faces() const;      ///< Return the faces of the StraightObject
            std::vector<Cell *> get_cells() const;      ///< Return the cells of the StraightObject

            size_t n_vertices() const; ///< Return the number of vertices of the StraightObject
            size_t n_edges() const;    ///< Return the number of edges of the StraightObject
            size_t n_faces() const;    ///< Return the number of faces of the StraightObject
            size_t n_cells() const;    ///< Return the number of cells of the StraightObject

            Vertex *vertex(const size_t i) const; ///< Return the i-th vertex of the StraightObject
            Edge *edge(const size_t i) const;     ///< Return the i-th edge of the StraightObject
            Face *face(const size_t i) const;     ///< Return the i-th face of the StraightObject
            Cell *cell(const size_t i) const;     ///< Return the i-th cell of the StraightObject

            void add_vertex(Vertex *vertex); ///< Add a vertex to the StraightObject
            void add_edge(Edge *edge);       ///< Add an edge to the StraightObject
            void add_face(Face *face);       ///< Add an face to the StraightObject
            void add_cell(Cell *cell);       ///< Add a cell to the StraightObject

            void remove_vertex(Vertex *vertex); ///< Remove a vertex from the StraightObject
            void remove_edge(Edge *edge);       ///< Remove a edge from the StraightObject
            void remove_face(Face *face);       ///< Remove a face from the StraightObject
            void remove_cell(Cell *cell);       ///< Remove a cell from the StraightObject

            int index_vertex(const Vertex *vertex) const; ///< Returns the local index of a vertex
            int index_edge(const Edge *edge) const;       ///< Returns the local index of an edge
            int index_face(const Face *face) const;       ///< Returns the local index of an face
            int index_cell(const Cell *cell) const;       ///< Returns the local index of a cell

            virtual bool test() const = 0; ///< Return a boolean determining if the geometries of the StraightObject are valid

        protected:
            size_t _index;
            VectorRd _center_mass;
            double _measure;
            double _diameter;
            bool _is_boundary;

            std::vector<Vertex *> _vertices;
            std::vector<Edge *> _edges;
            std::vector<Face *> _faces;
            std::vector<Cell *> _cells;
        };

        //@}
    } // end namespace StraightMesh

} // end namespace PolyMesh3D

#endif