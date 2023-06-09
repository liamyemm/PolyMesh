#include "QuadHandler.hpp"

namespace PolyMesh2D
{
    namespace Quadrature
    {
        template <typename MeshType>
        QuadHandler<MeshType>::QuadHandler(MeshType *mesh, size_t edge_doe, size_t cell_doe)
            : _edge_doe(edge_doe), _cell_doe(cell_doe)
        {
            // Reserve memory for edge_quads and cell_quads
            edge_quads.reserve(mesh->n_edges());
            cell_quads.reserve(mesh->n_cells());

            QuadratureRule<double> edge_legendre(gauss_jacobi(0, 0, _edge_doe));

            // Initialize edge quadrature rules
            for (auto &e : mesh->get_edges())
            {
                edge_quads.push_back(generate_quadrature_rule(e, edge_legendre));
            }

            QuadratureRule<double> cell_jacobi(gauss_jacobi(0, 1, _cell_doe));
            normalise_weights(cell_jacobi, 0, 1);

            // Initialize cell quadrature rules
            for (auto &c : mesh->get_cells())
            {
                // std::vector<QuadratureRule<double>> local_edge_quads(c->n_edges());
                std::vector<QuadratureRule<double>> local_edge_quads;
                local_edge_quads.reserve(c->n_edges());
                for (size_t iTF = 0; iTF < c->n_edges(); ++iTF)
                {
                    local_edge_quads.push_back(edge_quads[c->edge(iTF)->global_index()]);
                }
                cell_quads.push_back(generate_quadrature_rule(c, local_edge_quads, cell_jacobi));
            }
        }

        template <typename MeshType>
        double QuadHandler<MeshType>::integrate(const Functional::Function<2, 1> &func, typename MeshType::CellType *cell) const
        {
            typename MeshType::VertexType *vertex_split = nullptr;
            double beta = 0;
            if (func.poles.size() != 0)
            {
                // Can handle one singularity in the cell, located at a vertex.
                for (auto &pole : func.poles)
                {
                    for (auto &vert : cell->get_vertices())
                    {
                        if ((pole.location - vert->coords()).norm() < 1e-15)
                        {
                            vertex_split = vert;
                            beta = -(pole.order) + 1;
                            break;
                        }
                    }
                }
            }

            QuadratureRule<Eigen::Vector2d> cell_quad;
            if (vertex_split != nullptr)
            {
                QuadratureRule<double> cell_jacobi(gauss_jacobi(0, beta, _cell_doe));
                normalise_weights(cell_jacobi, 0, beta);

                std::vector<QuadratureRule<double>> local_edge_quads;
                local_edge_quads.reserve(cell->n_edges());
                for (size_t iTF = 0; iTF < cell->n_edges(); ++iTF)
                {
                    local_edge_quads.push_back(edge_quads[cell->edge(iTF)->global_index()]);
                }
                cell_quad = generate_quadrature_rule(cell, local_edge_quads, cell_jacobi, vertex_split);
            }
            else
            {
                cell_quad = cell_quads[cell->global_index()];
            }

            // Integrate function over cell
            double integral = 0.0;
            for (size_t i = 0; i < cell_quad.size(); i++)
            {
                Eigen::Vector2d point = cell_quad.point(i);
                double weight = cell_quad.weight(i);
                integral += weight * func.value(point);
            }
            return integral;
        }

        template <typename MeshType>
        double QuadHandler<MeshType>::integrate(const Functional::Function<1, 1> &func, typename MeshType::EdgeType *edge) const
        {
            bool at_tmin;
            double alf_or_bet = 0;
            if (func.poles.size() != 0)
            {
                // Can handle one singularity in the edge, located at a vertex.
                for (auto &pole : func.poles)
                {
                    if ((pole.location - edge->parameterisation().tmin) < 1e-15)
                    {
                        at_tmin = true;
                        alf_or_bet = -pole.order;
                        break;
                    }
                    else if ((pole.location - edge->parameterisation().tmax) < 1e-15)
                    {
                        at_tmin = false;
                        alf_or_bet = -pole.order;
                        break;
                    }
                }
            }

            QuadratureRule<double> edge_quad;

            if (alf_or_bet != 0)
            {
                QuadratureRule<double> edge_jacobi;
                if (at_tmin)
                {
                    edge_jacobi = gauss_jacobi(0, alf_or_bet, _edge_doe);
                    normalise_weights(edge_jacobi, 0, alf_or_bet);
                }
                else
                {
                    edge_jacobi = gauss_jacobi(alf_or_bet, 0, _edge_doe);
                    normalise_weights(edge_jacobi, alf_or_bet, 0);
                }
            }
            else
            {
                edge_quad = edge_quads[edge->global_index()];
            }

            // Integrate function over edge
            double integral = 0.0;
            for (size_t i = 0; i < edge_quad.size(); i++)
            {
                double point = edge_quad.point(i);
                double weight = edge_quad.weight(i);
                integral += weight * func.value(point);
            }
            return integral;
        }

        // Explicit instantiation of template class for supported mesh types
        // template class QuadHandler<StraightMesh::Mesh>;
        template class QuadHandler<PolyMesh2D::CurvedMesh::Mesh>;
    }
}