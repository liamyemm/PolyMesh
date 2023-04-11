#include "Enrichment.hpp"
#include <iostream>

namespace PolyMesh2D
{
    namespace HHOHELMHOLTZ
    {

        Enrichment::Enrichment(double lambda) : m_lambda(lambda) {}

        void Enrichment::globally_enrich(HybridCore &hho)
        {
            ScalarFunction2D func(this->glob_enrichment());
            ScalarFunction2D func_laplace(this->glob_enrichment_laplace());

            PolyMesh2D::CurvedMesh::Mesh *mesh = hho.get_mesh();

            for(size_t iT = 0; iT < mesh->n_cells(); ++iT)
            {
                hho.enrich_highorder_basis(iT, func);
                hho.enrich_cell_basis(iT, func_laplace);
            }

            for(size_t iE = 0; iE < mesh->n_edges(); ++iE)
            {
                ScalarFunction1D func_trace(this->glob_enrichment_neumann_trace(mesh->edge(iE)->parameterisation()));
                hho.enrich_edge_basis(iE, func_trace);
            }
        }

        void Enrichment::locally_enrich(HybridCore &hho)
        {

            PolyMesh2D::CurvedMesh::Mesh *mesh = hho.get_mesh();

            for(size_t iT = 0; iT < mesh->n_cells(); ++iT)
            {
                ScalarFunction2D func(this->loc_enrichment(mesh->cell(iT)->center_mass()));
                ScalarFunction2D func_laplace(this->loc_enrichment_laplace(mesh->cell(iT)->center_mass()));

                hho.enrich_highorder_basis(iT, func);
                hho.enrich_cell_basis(iT, func_laplace);
            }

            for(size_t iE = 0; iE < mesh->n_edges(); ++iE)
            {
                // if((iE != 0) && (iE != 1) && (iE != 2))
                // {
                //     continue;
                // }
                for(auto & cell : mesh->edge(iE)->get_cells())
                {
                    ScalarFunction1D func_trace(this->loc_enrichment_neumann_trace(cell->center_mass(), mesh->edge(iE)->parameterisation()));
                    hho.enrich_edge_basis(iE, func_trace);
                }
            }
        }

        ScalarFunction2D Enrichment::glob_enrichment() const
        {
            std::function<double(Eigen::Vector2d)> glob_enrich = [&](const Eigen::Vector2d x) -> double
            {
                return std::cos(m_lambda * x.norm());
            };

            std::function<Eigen::RowVector2d(Eigen::Vector2d)> D_glob_enrich = [&](const Eigen::Vector2d x) -> Eigen::RowVector2d
            {
                double multiplier = (x.norm() < 1E-15 ? 0.0 : 1.0 / x.norm());
                return -multiplier * m_lambda * std::sin(m_lambda * x.norm()) * x.transpose();
            };

            return ScalarFunction2D(glob_enrich, D_glob_enrich);
        }

        ScalarFunction2D Enrichment::glob_enrichment_laplace() const
        {
            std::function<double(Eigen::Vector2d)> glob_enrich_laplace = [&](const Eigen::Vector2d x) -> double
            {
                double r = x.norm();
                double val = -m_lambda * m_lambda * std::cos(m_lambda * r);
                if (r > 1E-15)
                {
                    val -= std::sin(m_lambda * r) / r;
                }
                return val;
            };

            std::function<Eigen::RowVector2d(Eigen::Vector2d)> D_glob_enrich_laplace = [&](const Eigen::Vector2d x) -> Eigen::RowVector2d
            {
                double r = x.norm();
                double k = m_lambda;
                if (r < 1E-15)
                {
                    return 0.0 * x.transpose();
                }
                return ((-k * r * std::cos(k * r) + (1.0 + std::pow(k, 3) * std::pow(r, 2)) * std::sin(k * r)) / std::pow(r, 3)) * x.transpose();
            };

            return ScalarFunction2D(glob_enrich_laplace, D_glob_enrich_laplace);
        }

        ScalarFunction1D Enrichment::glob_enrichment_neumann_trace(const Curve &edge_param) const
        {
            std::function<double(double)> neumann_trace = [edge_param, this](const double t) -> double
            {
                Eigen::Vector2d x = edge_param.value(t);
                if(x.norm() < 1E-15)
                {
                    return 0.0;
                }
                Eigen::Vector2d grad = -m_lambda * std::sin(m_lambda * x.norm()) * x / x.norm();
                return grad.dot(edge_param.normal(t));
            };

            return ScalarFunction1D(neumann_trace);
        }

        ScalarFunction2D Enrichment::loc_enrichment(const Eigen::Vector2d &cell_center) const
        {
            std::function<double(Eigen::Vector2d)> loc_enrich = [&](const Eigen::Vector2d x) -> double
            {
                Eigen::VectorXd y = x - cell_center;
                return std::cos(m_lambda * y.norm());
            };

            std::function<Eigen::RowVector2d(Eigen::Vector2d)> D_loc_enrich = [&](const Eigen::Vector2d x) -> Eigen::RowVector2d
            {
                Eigen::VectorXd y = x - cell_center;
                double multiplier = (y.norm() < 1E-15 ? 0.0 : 1.0 / y.norm());
                return -multiplier * m_lambda * std::sin(m_lambda * y.norm()) * y.transpose();
            };

            return ScalarFunction2D(loc_enrich, D_loc_enrich);
        }

        ScalarFunction2D Enrichment::loc_enrichment_laplace(const Eigen::Vector2d &cell_center) const
        {
            std::function<double(Eigen::Vector2d)> loc_enrich_laplace = [&](const Eigen::Vector2d x) -> double
            {
                Eigen::VectorXd y = x - cell_center;
                double r = y.norm();
                double val = -m_lambda * m_lambda * std::cos(m_lambda * r);
                if (r > 1E-15)
                {
                    val -= std::sin(m_lambda * r) / r;
                }
                return val;
            };

            std::function<Eigen::RowVector2d(Eigen::Vector2d)> D_loc_enrich_laplace = [&](const Eigen::Vector2d x) -> Eigen::RowVector2d
            {
                Eigen::VectorXd y = x - cell_center;
                double r = y.norm();
                double k = m_lambda;
                if (r < 1E-15)
                {
                    return 0.0 * x.transpose();
                }
                return ((-k * r * std::cos(k * r) + (1.0 + std::pow(k, 3) * std::pow(r, 2)) * std::sin(k * r)) / std::pow(r, 3)) * y.transpose();
            };
            return ScalarFunction2D(loc_enrich_laplace, D_loc_enrich_laplace);
        }

        ScalarFunction1D Enrichment::loc_enrichment_neumann_trace(const Eigen::Vector2d &cell_center, const Curve &edge_param) const
        {
            std::function<double(double)> neumann_trace = [edge_param, this, &cell_center](const double t) -> double
            {
                Eigen::Vector2d x = edge_param.value(t);
                Eigen::VectorXd y = x - cell_center;
                if(y.norm() < 1E-15)
                {
                    return 0.0;
                }
                Eigen::Vector2d grad = -m_lambda * std::sin(m_lambda * y.norm()) * y / y.norm();
                return grad.dot(edge_param.normal(t));
            };

            return ScalarFunction1D(neumann_trace);
        }

    }
}