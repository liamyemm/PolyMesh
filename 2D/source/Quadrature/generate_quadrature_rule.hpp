#include "QuadratureRule.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"
#include "Cell.hpp"

#ifndef _GENERATE_QUADRATURE_RULE_HPP
#define _GENERATE_QUADRATURE_RULE_HPP

namespace PolyMesh2D
{
    namespace Quadrature
    {
        /**
         * @brief Generates a quadrature rule for a curved edge using the provided 1D Gauss-Legendre quadrature rule.
         *
         * @param edge A pointer to the CurvedMesh::Edge the quadrature rule is generated for.
         * @param legendre_quad The 1D Gauss-Legendre quadrature rule to use.
         *
         * @return A quadrature rule consisting of nodes and weights.
         */
        const QuadratureRule<double> generate_quadrature_rule(const PolyMesh2D::CurvedMesh::Edge *edge, const QuadratureRule<double> &legendre_quad);

        // /**
        //  * @brief Generates a quadrature rule for a straight edge using the provided 1D Gauss-Legendre quadrature rule.
        //  *
        //  * @param edge A pointer to the StraightMesh::Edge the quadrature rule is generated for.
        //  * @param legendre_quad The 1D Gauss-Legendre quadrature rule to use.
        //  *
        //  * @return A quadrature rule consisting of nodes and weights.
        //  */
        // const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(const PolyMesh2D::StraightMesh::Edge *edge, const QuadratureRule<double> &legendre_quad);

        /**
         * @brief Generates a quadrature rule for a curved cell using the provided quadrature rules for its edges and center.
         *
         * @param cell A pointer to the CurvedMesh::Cell the quadrature rule is generated for.
         * @param edge_quads A vector of quadrature rules for the edges of the cell.
         * @param jacobi_quad The quadrature rule for the center of the cell.
         *
         * @return A quadrature rule consisting of nodes and weights.
         */
        const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(const PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads, const QuadratureRule<double> &jacobi_quad, const PolyMesh2D::CurvedMesh::Vertex *vertex_split = nullptr);

        // /**
        //  * @brief Generates a quadrature rule for a straight cell using the provided quadrature rules for its edges and center.
        //  *
        //  * @param cell A pointer to the StraightMesh::Cell the quadrature rule is generated for.
        //  * @param edge_quads A vector of quadrature rules for the edges of the cell.
        //  * @param jacobi_quad The quadrature rule for the center of the cell.
        //  *
        //  * @return A quadrature rule consisting of nodes and weights.
        //  */
        // const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(const PolyMesh2D::StraightMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads, const QuadratureRule<double> &jacobi_quad, const PolyMesh2D::CurvedMesh::Vertex *vertex_split = nullptr);

    }
}

#endif