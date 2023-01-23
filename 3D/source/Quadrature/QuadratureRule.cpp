#include "QuadratureRule.hpp"

namespace PolyMesh3D
{
    namespace Quadrature
    {
        const QuadratureRule generate_quadrature_rule(StraightMesh::Edge *edge, const GaussLegendre1D &quad)
        {
            QuadratureRule quad_rule_edge;

            for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
            {
                Eigen::Vector3d point = (edge->coords()[1] - edge->coords()[0]) * quad.point(iqn) + edge->coords()[0];
                double weight = edge->measure() * quad.weight(iqn);
                quad_rule_edge.push_back(QuadratureNode(point, weight));
            }

            return quad_rule_edge;
        }

        const QuadratureRule generate_quadrature_rule(StraightMesh::Face *face, const std::vector<QuadratureRule> &edge_quads, const GaussLegendre1D &quad)
        {
            QuadratureRule quad_rule_face;

            Eigen::Vector3d v0 = face->vertex(0)->coords();
            const size_t v0_index = face->vertex(0)->global_index();

            for (size_t iFE = 0; iFE < face->n_edges(); ++iFE)
            {
                StraightMesh::Edge *edge = face->edge(iFE);

                const size_t E_v0_index = edge->vertex(0)->global_index();
                const size_t E_v1_index = edge->vertex(1)->global_index();
                if( (E_v0_index == v0_index) || (E_v1_index == v0_index) )
                {
                    continue; // skip edges with vertices at x0
                }

                for (size_t iqn_edge = 0; iqn_edge < edge_quads[iFE].size(); ++iqn_edge)
                {
                    Eigen::Vector3d edge_point = edge_quads[iFE][iqn_edge].x;
                    double weighted_y_dot_nFE = edge_quads[iFE][iqn_edge].w * (edge_point - v0).dot(face->edge_orientation(iFE) * face->edge_normal(iFE));
                    for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
                    {
                        double t = quad.point(iqn);
                        Eigen::Vector3d face_point = t * edge_point + (1.0 - t) * v0;
                        double weight = t * quad.weight(iqn) * weighted_y_dot_nFE;
                        quad_rule_face.push_back(QuadratureNode(face_point, weight));
                    }
                }
            }

            return quad_rule_face;
        }

        const QuadratureRule generate_quadrature_rule(StraightMesh::Cell *cell, const std::vector<QuadratureRule> &face_quads, const GaussLegendre1D &quad)
        {
            QuadratureRule quad_rule_cell;

            Eigen::Vector3d v0 = cell->vertex(0)->coords();
            const size_t v0_index = cell->vertex(0)->global_index();

            for (size_t iTF = 0; iTF < cell->n_faces(); ++iTF)
            {
                StraightMesh::Face *face = cell->face(iTF);

                const size_t F_v0_index = face->vertex(0)->global_index();
                const size_t F_v1_index = face->vertex(1)->global_index();
                if( (F_v0_index == v0_index) || (F_v1_index == v0_index) )
                {
                    continue; // skip face with vertices at x0
                }

                for (size_t iqn_face = 0; iqn_face < face_quads[iTF].size(); ++iqn_face)
                {
                    Eigen::Vector3d face_point = face_quads[iTF][iqn_face].x;

                    double weighted_y_dot_nTF = face_quads[iTF][iqn_face].w * (face_point - v0).dot(cell->face_normal(iTF));
                    for (size_t iqn = 0; iqn < quad.n_points(); ++iqn)
                    {
                        double t = quad.point(iqn);
                        Eigen::Vector3d cell_point = t * face_point + (1.0 - t) * v0;
                        double weight = t * t * quad.weight(iqn) * weighted_y_dot_nTF;
                        quad_rule_cell.push_back(QuadratureNode(cell_point, weight));
                    }
                }
            }

            return quad_rule_cell;
        }
    }
}
