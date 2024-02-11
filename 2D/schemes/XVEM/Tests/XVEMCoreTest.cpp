#include "XVEMCoreTest.h"
#include "../PoissonSingularity.h"

namespace PolyMesh2D
{
    TEST_F(XVEMCoreTest, TestL2Projectors)
    {
        PolyMesh2D::XVEMCore xvem_l0_k2(curved_mesh.get(), 0, 2, true, true);

        std::function<double(double)> edge_const_lam = [](double x) -> double {return 3.3;};
        std::function<double(Eigen::Vector2d)> cell_const_lam = [](Eigen::Vector2d x) -> double {return 3.3;};

        PolyMesh2D::Functional::Function<1, 1> edge_const(edge_const_lam);
        PolyMesh2D::Functional::Function<2, 1> cell_const(cell_const_lam);

        for(size_t iE = 0; iE < curved_mesh->n_edges(); ++iE)
        {
            auto projection = xvem_l0_k2.l2_low_edge_projection(iE, edge_const);

            double val_at_center = 0.0;

            auto edge_param = curved_mesh->edge(iE)->parameterisation();
            double param_center = 0.5 * (edge_param.tmin + edge_param.tmax);
            for(size_t i = 0; i < xvem_l0_k2.local_low_edge_dofs(iE); ++i)
            {
                val_at_center += projection(i) * xvem_l0_k2.low_edge_basis(iE)->value(i, param_center);
            }

            ASSERT_NEAR(val_at_center, edge_const.value(param_center), 1E-13);
        }

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            auto projection = xvem_l0_k2.l2_low_cell_projection(iT, cell_const);

            double val_at_center = 0.0;

            for(size_t i = 0; i < xvem_l0_k2.local_low_cell_dofs(iT); ++i)
            {
                val_at_center += projection(i) * xvem_l0_k2.low_cell_basis(iT)->value(i, curved_mesh->cell(iT)->center_mass());
            }

            ASSERT_NEAR(val_at_center, cell_const.value(curved_mesh->cell(iT)->center_mass()), 1E-13);
        }

        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 2, 4, true, true);

        // PolyMesh2D::Functional::Function<2, 1> singularity = PolyMesh2D::PoissonSingularity(2.7);

        // for(size_t iT = 0; iT < xvem_l2_k4.get_mesh()->n_cells(); ++iT)
        // {
        //     xvem_l2_k4.enrich_cell_basis(iT, singularity);
        // }

        // for(size_t iE = 0; iE < xvem_l2_k4.get_mesh()->n_edges(); ++iE)
        // {
        //     xvem_l2_k4.enrich_edge_basis(iE, PolyMesh2D::Functional::trace(singularity, xvem_l2_k4.get_mesh()->edge(iE)->parameterisation()));
        // }

        std::function<double(double)> edge_quad_lam = [](double x) -> double {return 3.3 + 2.2 * x + x * x;};
        std::function<double(Eigen::Vector2d)> cell_quad_lam = [](Eigen::Vector2d x) -> double {return 3.3 + x(0) + 2.5 * x(1) + x(0) * x(1) + 3.0 * x(1) * x(1);};

        PolyMesh2D::Functional::Function<1, 1> edge_quad(edge_quad_lam);
        PolyMesh2D::Functional::Function<2, 1> cell_quad(cell_quad_lam);

        for(size_t iE = 0; iE < curved_mesh->n_edges(); ++iE)
        {
            auto projection = xvem_l2_k4.l2_low_edge_projection(iE, edge_quad);

            double val_at_t0 = 0.0;
            double val_at_t1 = 0.0;
            double val_at_center = 0.0;

            auto edge_param = curved_mesh->edge(iE)->parameterisation();
            double t0 = edge_param.tmin;
            double t1 = edge_param.tmax;
            double param_center = 0.5 * (t0 + t1);

            for(size_t i = 0; i < xvem_l2_k4.local_low_edge_dofs(iE); ++i)
            {
                val_at_t0 += projection(i) * xvem_l2_k4.low_edge_basis(iE)->value(i, t0);
                val_at_t1 += projection(i) * xvem_l2_k4.low_edge_basis(iE)->value(i, t1);
                val_at_center += projection(i) * xvem_l2_k4.low_edge_basis(iE)->value(i, param_center);
            }

            ASSERT_NEAR(val_at_t0, edge_quad.value(t0), 1E-12);
            ASSERT_NEAR(val_at_t1, edge_quad.value(t1), 1E-12);
            ASSERT_NEAR(val_at_center, edge_quad.value(param_center), 1E-12);
        }

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            auto projection = xvem_l2_k4.l2_low_cell_projection(iT, cell_quad);

            double val_at_vert0 = 0.0;
            double val_at_vert1 = 0.0;
            double val_at_vert2 = 0.0;
            double val_at_center = 0.0;

            Eigen::Vector2d vert0_coord = curved_mesh->cell(iT)->vertex(0)->coords();
            Eigen::Vector2d vert1_coord = curved_mesh->cell(iT)->vertex(1)->coords();
            Eigen::Vector2d vert2_coord = curved_mesh->cell(iT)->vertex(2)->coords();
            Eigen::Vector2d center_coord = curved_mesh->cell(iT)->center_mass();

            for(size_t i = 0; i < xvem_l2_k4.local_low_cell_dofs(iT); ++i)
            {
                val_at_vert0 += projection(i) * xvem_l2_k4.low_cell_basis(iT)->value(i, vert0_coord);
                val_at_vert1 += projection(i) * xvem_l2_k4.low_cell_basis(iT)->value(i, vert1_coord);
                val_at_vert2 += projection(i) * xvem_l2_k4.low_cell_basis(iT)->value(i, vert2_coord);
                val_at_center += projection(i) * xvem_l2_k4.low_cell_basis(iT)->value(i, center_coord);
            }

            ASSERT_NEAR(val_at_vert0, cell_quad.value(vert0_coord), 1E-12);
            ASSERT_NEAR(val_at_vert1, cell_quad.value(vert1_coord), 1E-12);
            ASSERT_NEAR(val_at_vert2, cell_quad.value(vert2_coord), 1E-12);
            ASSERT_NEAR(val_at_center, cell_quad.value(center_coord), 1E-12);
        }
    }

    TEST_F(XVEMCoreTest, TestEllipticProjector)
    {
        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 2, 3, true, false);

        PolyMesh2D::Functional::Function<2, 1> singularity = PolyMesh2D::PoissonSingularity(6.5);

        for(size_t iT = 0; iT < xvem_l2_k4.get_mesh()->n_cells(); ++iT)
        {
            xvem_l2_k4.enrich_cell_basis(iT, singularity);
        }

        for(size_t iE = 0; iE < xvem_l2_k4.get_mesh()->n_edges(); ++iE)
        {
            xvem_l2_k4.enrich_edge_basis(iE, PolyMesh2D::Functional::trace(singularity, xvem_l2_k4.get_mesh()->edge(iE)->parameterisation()));
        }

        std::function<double(Eigen::Vector2d)> cell_cubic_lam = [singularity](Eigen::Vector2d x) -> double {return singularity.value(x) + x(0) * x(1) + 3.0 * x(1) * x(1) + x(0) * x(0) * x(1) + 17.2 * x(1) * x(1) * x(0) - 12.7 * x(1) * x(1) * x(1) + 12.3 * x(0) * x(0) * x(0);};
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> cell_cubic_deriv = [singularity](Eigen::Vector2d x) -> Eigen::RowVector2d {return singularity.derivative(x) + Eigen::RowVector2d(x(1) + 2.0 * x(0) * x(1) + 17.2 * x(1) * x(1) + 3.0 * 12.3 * x(0) * x(0), x(0) + 6.0 * x(1) + x(0) * x(0) + 17.2 * 2.0 * x(1) * x(0) - 12.7 * 3.0 * x(1) * x(1));};

        // std::function<double(Eigen::Vector2d)> cell_cubic_lam = [](Eigen::Vector2d x) -> double {return x(0) * x(0) * x(0);};
        // std::function<Eigen::RowVector2d(Eigen::Vector2d)> cell_cubic_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(3.0 * x(0) * x(0), 0.0);};

        PolyMesh2D::Functional::Function<2, 1> cell_cubic(cell_cubic_lam, cell_cubic_deriv);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            auto projection = xvem_l2_k4.elliptic_projection(iT, cell_cubic);

            double val_at_vert0 = 0.0;
            double val_at_vert1 = 0.0;
            double val_at_vert2 = 0.0;
            double val_at_center = 0.0;

            Eigen::Vector2d vert0_coord = curved_mesh->cell(iT)->vertex(0)->coords();
            Eigen::Vector2d vert1_coord = curved_mesh->cell(iT)->vertex(1)->coords();
            Eigen::Vector2d vert2_coord = curved_mesh->cell(iT)->vertex(2)->coords();
            Eigen::Vector2d center_coord = curved_mesh->cell(iT)->center_mass();

            for(size_t i = 0; i < xvem_l2_k4.local_high_cell_dofs(iT); ++i)
            {
                val_at_vert0 += projection(i) * xvem_l2_k4.high_cell_basis(iT)->value(i, vert0_coord);
                val_at_vert1 += projection(i) * xvem_l2_k4.high_cell_basis(iT)->value(i, vert1_coord);
                val_at_vert2 += projection(i) * xvem_l2_k4.high_cell_basis(iT)->value(i, vert2_coord);
                val_at_center += projection(i) * xvem_l2_k4.high_cell_basis(iT)->value(i, center_coord);
            }

            ASSERT_NEAR(val_at_vert0, cell_cubic.value(vert0_coord), 1E-9);
            ASSERT_NEAR(val_at_vert1, cell_cubic.value(vert1_coord), 1E-9);
            ASSERT_NEAR(val_at_vert2, cell_cubic.value(vert2_coord), 1E-9);
            ASSERT_NEAR(val_at_center, cell_cubic.value(center_coord), 1E-9);
        }
    }

    TEST_F(XVEMCoreTest, TestInterpolate)
    {
        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 2, 4, false, true);

        std::function<double(Eigen::Vector2d)> sine_lam = [](Eigen::Vector2d x) -> double {return std::sin(x(0)) + std::cos(x(1));};
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> sine_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(std::cos(x(0)), -std::sin(x(1)));};

        PolyMesh2D::Functional::Function<2, 1> sine(sine_lam, sine_deriv);

        Eigen::VectorXd interp = xvem_l2_k4.interpolate(sine);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            Eigen::VectorXd cell_val = interp.segment(xvem_l2_k4.global_offset_T(iT), xvem_l2_k4.local_low_cell_dofs(iT));
            Eigen::VectorXd cell_proj = xvem_l2_k4.l2_low_cell_projection(iT, sine);
            ASSERT_NEAR(std::sqrt((cell_val - cell_proj).squaredNorm() / xvem_l2_k4.local_low_cell_dofs(iT)), 0.0, 1E-12);
        }

        for(size_t iE = 0; iE < curved_mesh->n_edges(); ++iE)
        {
            Eigen::VectorXd edge_val = interp.segment(xvem_l2_k4.global_offset_E(iE), xvem_l2_k4.local_low_edge_dofs(iE));
            Eigen::VectorXd edge_proj = xvem_l2_k4.l2_low_edge_projection(iE, PolyMesh2D::Functional::trace(sine, curved_mesh->edge(iE)->parameterisation()));
            ASSERT_NEAR(std::sqrt((edge_val - edge_proj).squaredNorm() / xvem_l2_k4.local_low_edge_dofs(iE)), 0.0, 1E-12);
        }

        for(size_t iV = 0; iV < curved_mesh->n_vertices(); ++iV)
        {
            double vertex_interp = interp(xvem_l2_k4.global_offset_V(iV));
            double vertex_val = sine.value(curved_mesh->vertex(iV)->coords());
            ASSERT_NEAR(vertex_interp, vertex_val, 1E-12);
        }
    }

    TEST_F(XVEMCoreTest, TestRestr)
    {
        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 2, 4, false, false);

        std::function<double(Eigen::Vector2d)> sine_lam = [](Eigen::Vector2d x) -> double {return std::sin(x(0)) + std::cos(x(1));};
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> sine_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(std::cos(x(0)), -std::sin(x(1)));};

        PolyMesh2D::Functional::Function<2, 1> sine(sine_lam, sine_deriv);

        Eigen::VectorXd interp = xvem_l2_k4.interpolate(sine);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            Eigen::VectorXd restriction = xvem_l2_k4.restr(interp, iT);

            Eigen::VectorXd cell_val = restriction.head(xvem_l2_k4.local_low_cell_dofs(iT));
            Eigen::VectorXd cell_proj = xvem_l2_k4.l2_low_cell_projection(iT, sine);
            ASSERT_NEAR(std::sqrt((cell_val - cell_proj).squaredNorm() / xvem_l2_k4.local_low_cell_dofs(iT)), 0.0, 1E-12);

            CurvedMesh::Cell *cell = curved_mesh->cell(iT);
            for(size_t iTE = 0; iTE < cell->n_edges(); ++iTE)
            {
                size_t iE = cell->edge(iTE)->global_index();
                Eigen::VectorXd edge_val = restriction.segment(xvem_l2_k4.local_offset_E(iT, iTE), xvem_l2_k4.local_low_edge_dofs(iE));
                Eigen::VectorXd edge_proj = xvem_l2_k4.l2_low_edge_projection(iE, PolyMesh2D::Functional::trace(sine, cell->edge(iTE)->parameterisation()));
                ASSERT_NEAR(std::sqrt((edge_val - edge_proj).squaredNorm() / xvem_l2_k4.local_low_edge_dofs(iE)), 0.0, 1E-12);
            }

            for(size_t iTV = 0; iTV < cell->n_vertices(); ++iTV)
            {
                double vertex_interp = restriction(xvem_l2_k4.local_offset_V(iT, iTV));
                double vertex_val = sine.value(cell->vertex(iTV)->coords());
                ASSERT_NEAR(vertex_interp, vertex_val, 1E-12);
            }
        }
    }

    TEST_F(XVEMCoreTest, TestEdgeRestrMat)
    {

        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 3, 3, true, true);

        PolyMesh2D::Functional::Function<2, 1> singularity = PolyMesh2D::PoissonSingularity(0.7);

        for(size_t iE = 0; iE < xvem_l2_k4.get_mesh()->n_edges(); ++iE)
        {
            xvem_l2_k4.enrich_edge_basis(iE, PolyMesh2D::Functional::trace(singularity, xvem_l2_k4.get_mesh()->edge(iE)->parameterisation()));
        }

        std::function<double(Eigen::Vector2d)> sine_lam = [](Eigen::Vector2d x) -> double {return std::sin(x(0)) + std::cos(x(1)) + 2.0;};
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> sine_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(std::cos(x(0)), -std::sin(x(1)));};

        PolyMesh2D::Functional::Function<2, 1> sine(sine_lam, sine_deriv);

        Eigen::VectorXd interp = xvem_l2_k4.interpolate(sine);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            Eigen::VectorXd restriction = xvem_l2_k4.restr(interp, iT);

            CurvedMesh::Cell *cell = curved_mesh->cell(iT);
            for(size_t iTE = 0; iTE < cell->n_edges(); ++iTE)
            {
                Eigen::VectorXd edge_restriction = xvem_l2_k4.edge_restrict_mat(iT, iTE) * restriction;

                size_t iE = cell->edge(iTE)->global_index();

                Eigen::VectorXd edge_val = edge_restriction.head(xvem_l2_k4.local_low_edge_dofs(iE));
                Eigen::VectorXd edge_proj = xvem_l2_k4.l2_low_edge_projection(iE, PolyMesh2D::Functional::trace(sine, cell->edge(iTE)->parameterisation()));
                ASSERT_NEAR(std::sqrt((edge_val - edge_proj).squaredNorm() / xvem_l2_k4.local_low_edge_dofs(iE)), 0.0, 1E-12);

                for(size_t iEV = 0; iEV < cell->edge(iTE)->n_vertices(); ++iEV)
                {
                    double vertex_interp = edge_restriction(xvem_l2_k4.local_low_edge_dofs(iE) + iEV);
                    double vertex_val = sine.value(cell->edge(iTE)->vertex(iEV)->coords());
                    ASSERT_NEAR(vertex_interp, vertex_val, 1E-12);
                }
            }
        }
    }

    TEST_F(XVEMCoreTest, TestContinuousReconstruction)
    {
        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 2, 4, true, true);

        PolyMesh2D::Functional::Function<2, 1> singularity = PolyMesh2D::PoissonSingularity(0.7);

        for(size_t iE = 0; iE < xvem_l2_k4.get_mesh()->n_edges(); ++iE)
        {
            xvem_l2_k4.enrich_edge_basis(iE, PolyMesh2D::Functional::trace(singularity, xvem_l2_k4.get_mesh()->edge(iE)->parameterisation()));
        }

        std::function<double(Eigen::Vector2d)> sine_lam = [](Eigen::Vector2d x) -> double {return std::sin(x(0)) + std::cos(x(1)) + 2.0;};
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> sine_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(std::cos(x(0)), -std::sin(x(1)));};

        PolyMesh2D::Functional::Function<2, 1> sine(sine_lam, sine_deriv);

        Eigen::VectorXd interp = xvem_l2_k4.interpolate(sine);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            Eigen::VectorXd restriction = xvem_l2_k4.restr(interp, iT);

            CurvedMesh::Cell *cell = curved_mesh->cell(iT);
            for(size_t iTE = 0; iTE < cell->n_edges(); ++iTE)
            {
                Eigen::VectorXd cont_recon = xvem_l2_k4.continuous_reconstruction(iT, iTE) * restriction;

                // project cont_recon onto low_edge_basis
                size_t iE = cell->edge(iTE)->global_index();

                Eigen::MatrixXd M_LOWE_MLOWE = xvem_l2_k4.get_quad_handle().l2_product(*xvem_l2_k4.low_edge_basis(iE), *xvem_l2_k4.low_edge_basis(iE), cell->edge(iTE), true);
                Eigen::MatrixXd M_LOWE_MHIGHE = xvem_l2_k4.get_quad_handle().l2_product(*xvem_l2_k4.low_edge_basis(iE), *xvem_l2_k4.high_edge_basis(iE), cell->edge(iTE), false);

                Eigen::VectorXd cont_recon_prog = M_LOWE_MLOWE.inverse() * M_LOWE_MHIGHE * cont_recon;

                Eigen::VectorXd edge_restriction = xvem_l2_k4.edge_restrict_mat(iT, iTE) * restriction;

                Eigen::VectorXd edge_val = edge_restriction.head(xvem_l2_k4.local_low_edge_dofs(iE));

                ASSERT_NEAR(std::sqrt((edge_val - cont_recon_prog).squaredNorm() / xvem_l2_k4.local_low_edge_dofs(iE)), 0.0, 1E-12);

                PolyMesh2D::Functional::Curve edge_param = cell->edge(iTE)->parameterisation();

                for(size_t iEV = 0; iEV < cell->edge(iTE)->n_vertices(); ++iEV)
                {
                    double vertex_interp = 0.0;

                    double t = edge_param.tmin;

                    if((edge_param.value(t) - cell->edge(iTE)->vertex(iEV)->coords()).norm() > 1E-14)
                    {
                        t = edge_param.tmax;
                    }

                    for(size_t i = 0; i < xvem_l2_k4.local_high_edge_dofs(iE); ++i)
                    {
                        vertex_interp += cont_recon(i) * xvem_l2_k4.high_edge_basis(iE)->value(i, t);
                    }
                    
                    
                    cont_recon(xvem_l2_k4.local_low_edge_dofs(iE) + iEV);
                    double vertex_val = sine_lam(curved_mesh->cell(iT)->edge(iTE)->vertex(iEV)->coords());
                    ASSERT_NEAR(vertex_interp, vertex_val, 1E-12);
                }
            }
        }
    }

    TEST_F(XVEMCoreTest, TestPotentialReconstruction)
    {
        PolyMesh2D::XVEMCore xvem_l2_k4(curved_mesh.get(), 0, 2, true, true);

        PolyMesh2D::Functional::Function<2, 1> singularity = PolyMesh2D::PoissonSingularity(0.7);

        for(size_t iT = 0; iT < xvem_l2_k4.get_mesh()->n_cells(); ++iT)
        {
            xvem_l2_k4.enrich_cell_basis(iT, singularity);
        }

        for(size_t iE = 0; iE < xvem_l2_k4.get_mesh()->n_edges(); ++iE)
        {
            xvem_l2_k4.enrich_edge_basis(iE, PolyMesh2D::Functional::trace(singularity, xvem_l2_k4.get_mesh()->edge(iE)->parameterisation()));
        }

        // std::function<double(Eigen::Vector2d)> cell_cubic_lam = [](Eigen::Vector2d x) -> double {return x(0) * x(1) + 3.0 * x(1) * x(1) + x(0) * x(0) * x(1) + 17.2 * x(1) * x(1) * x(0) - 12.7 * x(1) * x(1) * x(1) + 12.3 * x(0) * x(0) * x(0);};
        // std::function<Eigen::RowVector2d(Eigen::Vector2d)> cell_cubic_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(x(1) + 2.0 * x(0) * x(1) + 17.2 * x(1) * x(1) + 3.0 * 12.3 * x(0) * x(0), x(0) + 6.0 * x(1) + x(0) * x(0) + 17.2 * 2.0 * x(1) * x(0) - 12.7 * 3.0 * x(1) * x(1));};
        // }

        std::function<double(Eigen::Vector2d)> cell_cubic_lam = [](Eigen::Vector2d x) -> double {return x(0) * x(1) + 3.0 * x(1) * x(1) + x(0) * x(0);};
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> cell_cubic_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(x(1) + 2.0 * x(0), x(0) + 6.0 * x(1));};

        PolyMesh2D::Functional::Function<2, 1> cell_cubic(cell_cubic_lam, cell_cubic_deriv);
        Eigen::VectorXd interp = xvem_l2_k4.interpolate(cell_cubic);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            Eigen::VectorXd ellip = xvem_l2_k4.elliptic_projection(iT, cell_cubic);
            Eigen::VectorXd potential = xvem_l2_k4.potential_reconstruction(iT) * xvem_l2_k4.restr(interp, iT);

            ASSERT_NEAR(std::sqrt((ellip - potential).squaredNorm() / xvem_l2_k4.local_high_cell_dofs(iT)), 0.0, 1E-11);
        }

        interp = xvem_l2_k4.interpolate(singularity);

        for(size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            Eigen::VectorXd ellip = xvem_l2_k4.elliptic_projection(iT, singularity);
            Eigen::VectorXd potential = xvem_l2_k4.potential_reconstruction(iT) * xvem_l2_k4.restr(interp, iT);

            ASSERT_NEAR(std::sqrt((ellip - potential).squaredNorm() / xvem_l2_k4.local_high_cell_dofs(iT)), 0.0, 1E-10);
        }
    }
}

// Define the main function to run the tests
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}