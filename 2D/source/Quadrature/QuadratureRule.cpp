#include "QuadratureRule.hpp"
#include "Vertex.hpp"

namespace Quadrature
{
    const QuadratureRule<double> generate_quadrature_rule(PolyMesh2D::CurvedMesh::Edge *edge, const GaussLegendre1D &quad)
    {
        PolyMesh2D::Functional::Curve edge_param = edge->parameterisation();
        double a = edge_param.tmin;
        double b = edge_param.tmax;

        double len = b - a;
        assert(len > 0.0);

        QuadratureRule<double> quad_rule_edge;

        for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
        {
            double point = len * quad.point(iqn) + a;
            double weight = len * edge_param.derivative(point).norm() * quad.weight(iqn);
            quad_rule_edge.push_back(QuadratureNode<double>(point, weight));
            // std::cout << quad_rule_edge[iqn].x << "\n";
        }

        return quad_rule_edge;
    }

    const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads, const GaussLegendre1D &quad)
    {
        QuadratureRule<Eigen::Vector2d> quad_rule_cell;

        Eigen::Vector2d x0 = cell->vertex(0)->coords();
        // Eigen::Vector2d x0 = Eigen::Vector2d::Zero();

        const size_t v0_index = cell->vertex(0)->global_index();

        for (size_t iTF = 0; iTF < cell->n_edges(); ++iTF)
        {
            PolyMesh2D::CurvedMesh::Edge *edge = cell->edge(iTF);
             if(edge->is_straight())
             {
                 const size_t E_v0_index = edge->vertex(0)->global_index();
                 const size_t E_v1_index = edge->vertex(1)->global_index();
                 if( (E_v0_index == v0_index) || (E_v1_index == v0_index) )
                 {
                     continue; // skip straight edges with vertices at x0
                 }
             }
            for (size_t iqn_edge = 0; iqn_edge < edge_quads[iTF].size(); ++iqn_edge)
            {
                double edge_point = edge_quads[iTF][iqn_edge].x;
                // double edge_weight = edge_quads[iTF][iqn_edge].w;

                Eigen::Vector2d x = edge->parameterisation().value(edge_point);
                double weighted_y_dot_nTF = edge_quads[iTF][iqn_edge].w * (x - x0).dot(cell->edge_normal(iTF, edge_point));
                for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
                {
                    double t = quad.point(iqn);
                    Eigen::Vector2d cell_point = t * x + (1.0 - t) * x0;
                    double weight = t * quad.weight(iqn) * weighted_y_dot_nTF;
                    quad_rule_cell.push_back(QuadratureNode<Eigen::Vector2d>(cell_point, weight));
                }
            }
        }

        return quad_rule_cell;
    }

    const QuadratureRule<Eigen::Vector2d> generate_weighted_boundary_rule(PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads)
    {
        QuadratureRule<Eigen::Vector2d> quad_rule_cell;

        for (size_t iTF = 0; iTF < cell->n_edges(); ++iTF)
        {
            for (size_t iqn_edge = 0; iqn_edge < edge_quads[iTF].size(); ++iqn_edge)
            {
                double edge_point = edge_quads[iTF][iqn_edge].x;
                Eigen::Vector2d point = cell->edge(iTF)->parameterisation().value(edge_point);
                double weight = edge_quads[iTF][iqn_edge].w * point.dot(cell->edge_normal(iTF, edge_point));
                quad_rule_cell.push_back(QuadratureNode<Eigen::Vector2d>(point, weight));
            }
        }

        return quad_rule_cell;
    }
}
