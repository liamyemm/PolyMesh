#include "generate_quadrature_rule.hpp"

namespace PolyMesh2D
{
    namespace Quadrature
    {
        const QuadratureRule<double> generate_quadrature_rule(const PolyMesh2D::CurvedMesh::Edge *edge, const QuadratureRule<double> &quad)
        {
            PolyMesh2D::Functional::Curve edge_param = edge->parameterisation();
            double a = edge_param.tmin;
            double b = edge_param.tmax;

            double len = b - a;
            assert(len > 0.0);

            QuadratureRule<double> quad_rule_edge;

            std::vector<double> points;
            std::vector<double> weights;

            for (size_t iqn = 0; iqn < quad.size(); ++iqn)
            {
                double point = len * quad.point(iqn) + a;
                double weight = len * edge_param.derivative(point).norm() * quad.weight(iqn);
                points.push_back(point);
                weights.push_back(weight);
            }

            quad_rule_edge.points = points;
            quad_rule_edge.weights = weights;

            return quad_rule_edge;
        }

        // const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(PolyMesh2D::StraightMesh::Edge *edge, const QuadratureRule<double> &quad)
        // {
        //     Eigen::Vector2d x0 = edge->vertex(0)->coords();
        //     Eigen::Vector2d x1 = edge->vertex(1)->coords();

        //     QuadratureRule<Eigen::Vector2d> quad_rule_edge;

        //     std::vector<Eigen::Vector2d> points;
        //     std::vector<double> weights;

        //     for (size_t iqn = 0; iqn < quad.size(); ++iqn)
        //     {
        //         Eigen::Vector2d point = (x1 - x0) * quad.point(iqn) + x0;
        //         double weight = (x1 - x0).norm() * quad.weight(iqn);
        //         points.push_back(point);
        //         weights.push_back(weight);
        //     }

        //     quad_rule_edge.points = points;
        //     quad_rule_edge.weights = weights;

        //     return quad_rule_edge;
        // }


        const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(const PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads, const QuadratureRule<double> &quad, const PolyMesh2D::CurvedMesh::Vertex *vertex_split)
        {
            Eigen::Vector2d x0 = cell->vertex(0)->coords();
            size_t v0_index = cell->vertex(0)->global_index();

            if(vertex_split != nullptr)
            {
                std::vector<PolyMesh2D::CurvedMesh::Vertex *> cell_verts = cell->get_vertices();
                assert(std::find(cell_verts.begin(), cell_verts.end(), vertex_split) != cell_verts.end());

                x0 = vertex_split->coords();
                v0_index = vertex_split->global_index();
            }

            QuadratureRule<Eigen::Vector2d> quad_rule_cell;

            std::vector<Eigen::Vector2d> points;
            std::vector<double> weights;

            for (size_t iTF = 0; iTF < cell->n_edges(); ++iTF)
            {
                PolyMesh2D::CurvedMesh::Edge *edge = cell->edge(iTF);
                if (edge->is_straight())
                {
                    const size_t E_v0_index = edge->vertex(0)->global_index();
                    const size_t E_v1_index = edge->vertex(1)->global_index();
                    if ((E_v0_index == v0_index) || (E_v1_index == v0_index))
                    {
                        continue; // skip straight edges with vertices at x0
                    }
                }
                for (size_t iqn_edge = 0; iqn_edge < edge_quads[iTF].size(); ++iqn_edge)
                {
                    double edge_point = edge_quads[iTF][iqn_edge].x;

                    Eigen::Vector2d x = edge->parameterisation().value(edge_point);
                    double weighted_y_dot_nTF = edge_quads[iTF][iqn_edge].w * (x - x0).dot(cell->edge_normal(iTF, edge_point));
                    for (size_t iqn = 0; iqn < quad.size(); ++iqn)
                    {
                        double t = quad.point(iqn);
                        Eigen::Vector2d cell_point = t * x + (1.0 - t) * x0;
                        double weight = t * quad.weight(iqn) * weighted_y_dot_nTF;
                        points.push_back(cell_point);
                        weights.push_back(weight);
                    }
                }
            }

            quad_rule_cell.points = points;
            quad_rule_cell.weights = weights;

            return quad_rule_cell;
        }

        // const QuadratureRule<Eigen::Vector2d> generate_quadrature_rule(PolyMesh2D::StraightMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads, const QuadratureRule<double> &quad, const PolyMesh2D::StraightMesh::Vertex *vertex_split)
        // {

        //     Eigen::Vector2d x0 = cell->vertex(0)->coords();
        //     const size_t v0_index = cell->vertex(0)->global_index();

        //     if(vertex_split != nullptr)
        //     {
        //         std::vertex<PolyMesh2D::StraightMesh::Vertex *> cell_verts = cell->get_vertices();
        //         assert(std::find(cell_verts.begin(), cell_verts.end(), vertex_split) != cell_verts.end());

        //         x0 = vertex_split->coords();
        //         v0_index = vertex_split->global_index();
        //     }

        //     QuadratureRule<Eigen::Vector2d> quad_rule_cell;

        //     std::vector<Eigen::Vector2d> points;
        //     std::vector<double> weights;

        //     for (size_t iTF = 0; iTF < cell->n_edges(); ++iTF)
        //     {
        //         PolyMesh2D::StraightMesh::Edge *edge = cell->edge(iTF);
        //         const size_t E_v0_index = edge->vertex(0)->global_index();
        //         const size_t E_v1_index = edge->vertex(1)->global_index();
        //         if ((E_v0_index == v0_index) || (E_v1_index == v0_index))
        //         {
        //             continue; // skip straight edges with vertices at x0
        //         }
        //         for (size_t iqn_edge = 0; iqn_edge < edge_quads[iTF].size(); ++iqn_edge)
        //         {
        //             Eigen::Vector2d x = edge_quads[iTF][iqn_edge].x;
        //             double weighted_y_dot_nTF = edge_quads[iTF][iqn_edge].w * (x - x0).dot(cell->edge_normal(iTF));
        //             for (size_t iqn = 0; iqn < quad.size(); ++iqn)
        //             {
        //                 double t = quad.point(iqn);
        //                 Eigen::Vector2d cell_point = t * x + (1.0 - t) * x0;
        //                 double weight = t * quad.weight(iqn) * weighted_y_dot_nTF;
        //                 points.push_back(cell_point);
        //                 weights.push_back(weight);
        //             }
        //         }
        //     }

        //     quad_rule_cell.points = points;
        //     quad_rule_cell.weights = weights;

        //     return quad_rule_cell;
        // }

        // const QuadratureRule<Eigen::Vector2d> generate_weighted_boundary_rule(PolyMesh2D::CurvedMesh::Cell *cell, const std::vector<QuadratureRule<double>> &edge_quads)
        // {
        //     QuadratureRule<Eigen::Vector2d> quad_rule_cell;

        //     for (size_t iTF = 0; iTF < cell->n_edges(); ++iTF)
        //     {
        //         for (size_t iqn_edge = 0; iqn_edge < edge_quads[iTF].size(); ++iqn_edge)
        //         {
        //             double edge_point = edge_quads[iTF][iqn_edge].x;
        //             Eigen::Vector2d point = cell->edge(iTF)->parameterisation().value(edge_point);
        //             double weight = edge_quads[iTF][iqn_edge].w * point.dot(cell->edge_normal(iTF, edge_point));
        //             quad_rule_cell.push_back(QuadratureNode<Eigen::Vector2d>(point, weight));
        //         }
        //     }

        //     return quad_rule_cell;
        // }
    }
}
