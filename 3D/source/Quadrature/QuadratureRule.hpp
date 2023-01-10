#include "GaussLegendre.hpp"
#include "Edge.hpp"
#include "Cell.hpp"

#ifndef _QUADRATURERULE_HPP
#define _QUADRATURERULE_HPP

namespace Quadrature
{
    template <typename point_type> // doubles or Eigen::Vector2d
    struct QuadratureNode
    {
        point_type x;
        double w;
        QuadratureNode(point_type x, double w) : x(x), w(w) {}
    };

    template <typename point_type>
    using QuadratureRule = std::vector<QuadratureNode<point_type>>;

    const QuadratureRule<double> generate_quadrature_rule(PolyMesh2D::CurvedMesh::Edge *edge, const GaussLegendre1D &quad);
    const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads, const GaussLegendre1D &quad);
    const QuadratureRule<Eigen::Vector2d> generate_weighted_boundary_rule(PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads);
}
#endif