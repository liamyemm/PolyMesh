#include <Polytope.hpp>

#ifndef _CELL_HPP
#define _CELL_HPP

namespace PolyMesh2D
{
    namespace CurvedMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Cell : public Polytope
        {
        public:
            Cell(size_t index, std::vector<Edge *> edges); ///< Constructor for a Cell. Calls the base class constructor, sets _edges to @arg edges and computes the cell data (center_mass, diameter, measure, edge orientations)

            int edge_orientation(const size_t edge_index) const;           ///< Return the orientation of the Edge located at @arg edge_index.
            VectorRd edge_normal(const size_t edge_index, double t) const; ///< Return the outer normal of the Edge located at @arg edge_index evaluated at the parameter @t.

            bool is_straight() const; ///< Return a boolean stating whether or not all the edges of the Cell are straight.

            bool test() const;

            void set_edge_orientations();

        private:
            bool _is_straight;
            std::vector<int> _edge_orientations;
        };

        //@}
    }

    namespace StraightMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        using Triangle = std::array<VectorRd, 3>;

        double signed_area(VectorRd A, VectorRd B, VectorRd C);

        class Cell : public Polytope
        {
        public:
            Cell(size_t index, std::vector<Triangle> triangles); ///< Constructor for a Cell. Calls the base class constructor, sets _triangles to @arg triangles and computes the cell data (center_mass, diameter, measure, edge orientations)

            std::vector<Triangle> triangulation() const { return _triangles; } ///< Return the triangulation of the Polytope

            int edge_orientation(const size_t edge_index) const;           ///< Return the orientation of the Edge located at @arg edge_index.
            VectorRd edge_normal(const size_t edge_index) const; ///< Return the outer normal of the Edge located at @arg edge_index.

            void plot_triangulation(std::ofstream *out) const; ///< Plot the triangles to @arg out

        private:
            std::vector<Triangle> _triangles;
            std::vector<int> _edge_orientations;

            void set_edge_orientations();
        };

        //@}
    }
}

#endif