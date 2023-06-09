#include <StokesCore.hpp>
#include <parallel_for.hpp>
#include "common.hpp"
#include "vtu_writer.hpp"
#include <iomanip>

using namespace PolyMesh2D;

Eigen::MatrixXd orthonormalize(Eigen::MatrixXd &gram_mat)
{
    // size_t dimension = gram_mat.rows();
    // assert((size_t)gram_mat.cols() == dimension);

    // Eigen::MatrixXd B = Eigen::MatrixXd::Zero(dimension, dimension);
    // double norm = std::sqrt(gram_mat(0, 0));
    // B(0, 0) = (1.0 / norm);
    // gram_mat.row(0) /= norm;
    // Eigen::VectorXd tmp = gram_mat.row(0).transpose();
    // gram_mat.col(0) = tmp;

    // for (size_t i = 1; i < dimension; ++i)
    // {
    //     Eigen::RowVectorXd coeffs = Eigen::RowVectorXd::Zero(i);
    //     for (size_t j = 0; j < i; j++)
    //     {
    //         coeffs(j) = -gram_mat(i, j);
    //     }

    //     for (size_t j = 0; j < i; j++)
    //     {
    //         gram_mat.row(i) += coeffs(j) * gram_mat.row(j);
    //     }

    //     norm = std::sqrt(gram_mat(i, i));
    //     coeffs /= norm;
    //     gram_mat.row(i) /= norm;

    //     tmp = gram_mat.row(i).transpose();
    //     gram_mat.col(i) = tmp;

    //     B.row(i).head(i) = coeffs * B.topLeftCorner(i, i);
    //     B(i, i) = (1.0 / norm);
    // }
    // return B;

    Eigen::LLT<Eigen::MatrixXd> llt(gram_mat);
    Eigen::MatrixXd L = llt.matrixL();
    return L.inverse();
}

StokesCore::StokesCore(
    MeshType *mesh_ptr, ///< A pointer to the loaded mesh
    size_t cell_deg,    ///< The degree of the cell polynomials
    size_t edge_deg,    ///< The degree of the edge polynomials
    bool use_threads,   ///< Optional argument to indicate if threads should be used
    bool ortho)

    : m_mesh(mesh_ptr),
      m_cell_deg(cell_deg),
      m_edge_deg(edge_deg),
      // m_quad_handle(m_mesh, 2 * m_edge_deg + 2, 2 * m_edge_deg + 2),
      m_quad_handle(m_mesh, 45, 45),
      m_use_threads(use_threads),
      m_ortho(ortho),
      m_highorder_basis(mesh_ptr->n_cells()),
      m_cell_basis(mesh_ptr->n_cells()),
      m_edge_basis(mesh_ptr->n_edges()),
      m_pressure_basis(mesh_ptr->n_cells())
{
    // boost::timer::cpu_timer construction_timer;
    // construction_timer.start();

    // std::cout << "\n[StokesCore] Constructing bases\n";

    // Create extended bases
    std::function<void(size_t, size_t)> construct_all_highorder_bases = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            m_highorder_basis[iT].reset(new CellBasisType(construct_cell_basis(iT, m_edge_deg + 1)));
        }
    };
    parallel_for(m_mesh->n_cells(), construct_all_highorder_bases, m_use_threads);

    // Create cell bases
    std::function<void(size_t, size_t)> construct_all_cell_bases = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            m_cell_basis[iT].reset(new CellBasisType(construct_cell_basis(iT, m_cell_deg)));
        }
    };
    parallel_for(m_mesh->n_cells(), construct_all_cell_bases, m_use_threads);

    // Create edge bases
    std::function<void(size_t, size_t)> construct_all_edge_bases = [&](size_t start, size_t end) -> void
    {
        for (size_t iF = start; iF < end; iF++)
        {
            m_edge_basis[iF].reset(new EdgeBasisType(construct_edge_basis(iF)));
        }
    };
    parallel_for(m_mesh->n_edges(), construct_all_edge_bases, m_use_threads);

    // Create pressure bases
    std::function<void(size_t, size_t)> construct_all_pressure_bases = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            m_pressure_basis[iT].reset(new ScalarBasisType(construct_pressure_basis(iT, m_edge_deg)));
        }
    };
    parallel_for(m_mesh->n_cells(), construct_all_pressure_bases, m_use_threads);

    // construction_timer.stop();
    // std::cout << "     Construction time = " << construction_timer.elapsed().wall * 1E-9 << "s\n";
}

StokesCore::ScalarBasisType StokesCore::construct_pressure_basis(const size_t iT, const size_t degree) const
{
    Eigen::Vector2d xT = m_mesh->cell(iT)->center_mass();
    double hT = m_mesh->cell(iT)->diam();

    std::function<Eigen::Vector2d(Eigen::Vector2d)> transform_val = [xT, hT](Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return (1.0 / hT) * (x - xT);
    };

    std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [xT, hT](Eigen::Vector2d x) -> Eigen::Matrix2d
    {
        return (1.0 / hT) * Eigen::Matrix2d::Identity();
    };

    Functional::VectorFunction2D transform(transform_val, transform_deriv);

    // Basis of monomials
    Functional::ScalarBasis2D mono_basis = Functional::MonomialScalarBasis(transform, degree);

    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->cell(iT), true);

    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension());

    if (m_ortho)
    {
        B = orthonormalize(basis_integrals);
    }

    return ScalarBasisType(mono_basis, B);
}

StokesCore::CellBasisType StokesCore::construct_cell_basis(const size_t iT, const size_t degree) const
{
    Eigen::Vector2d xT = m_mesh->cell(iT)->center_mass();
    double hT = m_mesh->cell(iT)->diam();

    std::function<Eigen::Vector2d(Eigen::Vector2d)> transform_val = [xT, hT](Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return (1.0 / hT) * (x - xT);
    };

    std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [xT, hT](Eigen::Vector2d x) -> Eigen::Matrix2d
    {
        return (1.0 / hT) * Eigen::Matrix2d::Identity();
    };

    Functional::VectorFunction2D transform(transform_val, transform_deriv);

    // Basis of monomials
    Functional::VectorBasis2D mono_basis = Functional::MonomialVectorBasis(transform, degree);

    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->cell(iT), true);

    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension());

    if (m_ortho)
    {
        B = orthonormalize(basis_integrals);
    }

    return CellBasisType(mono_basis, B);
}

void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}

void removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

    if (colToRemove < numCols)
        matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

    matrix.conservativeResize(numRows, numCols);
}

StokesCore::EdgeBasisType StokesCore::construct_edge_basis(const size_t iE) const
{
    PolyMesh2D::Functional::VectorBasis1D mono_basis;
    if (m_mesh->edge(iE)->is_straight())
    {
        mono_basis = PolyMesh2D::Functional::MonomialVectorBasis(m_mesh->edge(iE)->parameterisation(), m_edge_deg);
    }
    else
    {
        mono_basis = PolyMesh2D::Functional::CurvedEdgeVectorBasis(m_mesh->edge(iE)->parameterisation(), m_edge_deg);
    }

    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->edge(iE), true);

    Eigen::FullPivLU<Eigen::MatrixXd> lu(basis_integrals);

    double threshold = 1E-15;
    lu.setThreshold(threshold);

    Eigen::MatrixXd image = lu.image(basis_integrals);

    Eigen::MatrixXd new_basis_integrals = basis_integrals;

    if (image.cols() != basis_integrals.cols())
    {
        // std::cout << "\n\n" << lu.rank() << "\n";
        // std::cout << mono_basis.dimension() << "\n";
        int shift = 0;
        for (int basis_col = 0; basis_col < basis_integrals.cols(); ++basis_col)
        {
            bool to_remove = true;
            for (int image_col = 0; image_col < image.cols(); ++image_col)
            {
                // if ((image.col(image_col) - basis_integrals.col(basis_col)).norm() < threshold)
                // {
                if (image.col(image_col) == basis_integrals.col(basis_col))
                {
                    removeColumn(image, image_col);
                    to_remove = false;
                    break;
                }
            }
            if (to_remove)
            {
                mono_basis.remove_basis_function(basis_col - shift);
                ++shift;
            }
        }
        // std::cout << mono_basis.dimension() << "\n";

        new_basis_integrals = m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->edge(iE), true);
    }
    else
    {
        new_basis_integrals = basis_integrals;
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension());

    if (m_ortho)
    {
        B = orthonormalize(new_basis_integrals);
    }

    return EdgeBasisType(mono_basis, B);
}

PolyMesh2D::Functional::ScalarFamily1D StokesCore::construct_scalar_edge_basis(const size_t iE, const size_t degree) const
{
    assert(m_mesh->edge(iE)->is_straight());

    PolyMesh2D::Functional::ScalarBasis1D mono_basis;
    mono_basis = PolyMesh2D::Functional::MonomialScalarBasis(m_mesh->edge(iE)->parameterisation(), degree);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->edge(iE), true);

    Eigen::MatrixXd B = orthonormalize(basis_integrals);

    return PolyMesh2D::Functional::ScalarFamily1D(mono_basis, B);
}

Functional::Function<1, 1> ramp()
{
    std::function<double(double)> ramp = [](double x) -> double
    {
        if (x < 0.027) // difference between function and zero is numerically zero below this threshold
        {
            return 0.0;
        }
        if (x > 0.973) // difference between function and one is numerically zero above this threshold
        {
            return 1.0;
        }
        return 1.0 / (1.0 + std::exp(1.0 / x - 1.0 / (1.0 - x)));
    };
    std::function<double(double)> ramp_deriv = [](double x) -> double
    {
        if (x < 0.022) // difference between function and zero is numerically zero below this threshold
        {
            return 0.0;
        }
        if (x > 0.978) // difference between function and zero is numerically zero above this threshold
        {
            return 0.0;
        }
        return (std::exp(1.0 / (x - x * x)) * (1.0 - 2.0 * (1 - x) * x)) / (std::pow((std::exp(1.0 / (1.0 - x)) + std::exp(1.0 / x)) * x * (1 - x), 2));
    };

    return Functional::Function<1, 1>(ramp, ramp_deriv);
}

Functional::Function<1, 1> ramp_shift(double eps1, double eps2)
{
    assert(eps1 >= 0.0 && eps1 < 1.0);

    assert(eps2 > 0.0 && eps2 <= 1.0);

    // eps1 determines right shift - point at which function starts ramping up
    // eps2 determines scale over which ramp up occurs

    std::function<double(double)> shift = [eps1, eps2](double x) -> double
    {
        return (x - eps1) / eps2;
    };
    std::function<double(double)> shift_deriv = [eps2](double x) -> double
    {
        return x / eps2;
    };

    return Functional::Function<1, 1>(shift, shift_deriv);
}

Functional::Function<2, 1> reference_bubble(double eps1, double eps2)
{
    auto ramp_func = Functional::compose(ramp(), ramp_shift(eps1, eps2));

    std::function<double(Eigen::Vector2d)> bubble = [ramp_func](Eigen::Vector2d x) -> double
    {
        return ramp_func.value(x(0)) * ramp_func.value(x(1)) * ramp_func.value(1.0 - x(0) - x(1));
    };
    std::function<Eigen::RowVector2d(Eigen::Vector2d)> bubble_deriv = [ramp_func](Eigen::Vector2d x) -> Eigen::RowVector2d
    {
        double d1 = ramp_func.value(x(1)) * (ramp_func.value(1.0 - x(0) - x(1)) * ramp_func.derivative(x(0)) - ramp_func.value(x(0)) * ramp_func.derivative(1.0 - x(0) - x(1)));
        double d2 = ramp_func.value(x(0)) * (ramp_func.value(1.0 - x(0) - x(1)) * ramp_func.derivative(x(1)) - ramp_func.value(x(1)) * ramp_func.derivative(1.0 - x(0) - x(1)));
        return Eigen::RowVector2d(d1, d2);
    };
    return Functional::Function<2, 1>(bubble, bubble_deriv);
}

Functional::Function<2, 1> bubble(const CurvedMesh::Cell *cell, double eps1, double eps2)
{
    assert(cell->n_vertices() == 3);

    auto ref_bubble = reference_bubble(eps1, eps2);

    Eigen::Matrix2d A;
    A.col(0) = cell->vertex(1)->coords() - cell->vertex(0)->coords();
    A.col(1) = cell->vertex(2)->coords() - cell->vertex(0)->coords();

    Eigen::Matrix2d B = A.inverse();
    Eigen::Vector2d b = cell->vertex(0)->coords();

    std::function<Eigen::Vector2d(Eigen::Vector2d)> transform = [B, b](Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return B * (x - b);
    };
    std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [B](Eigen::Vector2d x) -> Eigen::Matrix2d
    {
        return B;
    };
    return Functional::compose(ref_bubble, Functional::Function<2, 2>(transform, transform_deriv));
}

Eigen::VectorXd StokesCore::pressure_robust_RHS(const size_t iT, const Functional::VectorFunction2D &source, const std::vector<Functional::Function<2, 1>> &RTN_enrichment, const std::vector<Functional::Function<2, 2>> &laplace_cell_enrichment) const
{
    assert(m_mesh->cell(iT)->n_edges() == 3);
    // generate an RTN basis
    Eigen::Vector2d xT = m_mesh->cell(iT)->center_mass();
    double hT = m_mesh->cell(iT)->diam();

    assert(RTN_enrichment.size() == laplace_cell_enrichment.size());

    // auto bubble_enrich = Functional::curl(the_bubble * func);

    // auto func = RTN_enrichment[0];
    // auto the_bubble = reference_bubble(0, 0.3);
    // // auto bubble_enrich = the_bubble * func;

    // Eigen::VectorXd xcoord = Eigen::VectorXd::Zero(m_mesh->n_vertices());
    // Eigen::VectorXd ycoord = Eigen::VectorXd::Zero(m_mesh->n_vertices());

    // for (size_t iV = 0; iV < m_mesh->n_vertices(); iV++)
    // {
    //     xcoord(iV) = bubble_enrich.value(m_mesh->vertex(iV)->coords())(0);

    //     // ycoord(iV) = bubble_enrich.value(m_mesh->vertex(iV)->coords())(1);
    // }

    // VtuWriter plotdata(m_mesh);

    // // plotdata.write_to_vtu("exact-" + plot_file + ".vtu", exact);
    // plotdata.write_to_vtu("bubblex.vtu", xcoord);
    // // plotdata.write_to_vtu("bubbley.vtu", ycoord);

    // exit(1);

    std::function<Eigen::Vector2d(Eigen::Vector2d)> transform_val = [xT, hT](Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return (1.0 / hT) * (x - xT);
    };

    std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [xT, hT](Eigen::Vector2d x) -> Eigen::Matrix2d
    {
        return (1.0 / hT) * Eigen::Matrix2d::Identity();
    };

    Functional::VectorFunction2D transform(transform_val, transform_deriv);

    Functional::VectorBasis2D RTN_basis_non_ortho(RaviartThomasNedelecBasis(transform, m_edge_deg + 1));

    // std::ofstream bubble_out("bubble_plot.tsv");

    for (auto &func : RTN_enrichment)
    {
        // auto the_bubble = reference_bubble(0, 0.2);
        auto the_bubble = bubble(m_mesh->cell(iT), 0.1, 0.0002);
        auto bubble_enrich = Functional::curl(the_bubble * func);

        std::vector<Functional::Pole<Eigen::Vector2d>> poles;
        bubble_enrich.set_poles(poles); // empty poles
        // auto bubble_enrich = Functional::curl(the_bubble);

        // std::cout << "\n" << bubble_enrich.value(m_mesh->cell(iT)->center_mass()) << "\n";

        RTN_basis_non_ortho.add_basis_function(bubble_enrich);
        // for (size_t jT = 0; jT < m_mesh->n_cells(); ++jT)
        // {
        //     for (auto & point : m_quad_handle.get_cell_quad(jT).points)
        //     {
        //         bubble_out << std::setprecision(12) << std::setw(20) << std::left << point(0) << std::setprecision(12) << std::setw(20) << std::left << point(1) << std::setprecision(12) << std::setw(20) << the_bubble.value(point) << std::endl;
        //     }
        // }
    }
    // bubble_out.close();

    // exit(1);

    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(RTN_basis_non_ortho, RTN_basis_non_ortho, m_mesh->cell(iT));
    Eigen::MatrixXd B = orthonormalize(basis_integrals);

    CellBasisType RTN_basis(RTN_basis_non_ortho, B);

    // Eigen::MatrixXd RTN_mat = Eigen::MatrixXd::Zero(RTN_basis.dimension(), local_cell_dofs(iT) + local_boundary_dofs(iT));
    // Eigen::MatrixXd RTN_mass_mat = Eigen::MatrixXd::Zero(cell_test_space.dimension() + bdry_test_space.dimension(), RTN_basis.dimension());
    // Eigen::MatrixXd RTN_RHS_mat = Eigen::MatrixXd::Zero(cell_test_space.dimension() + bdry_test_space.dimension(), local_cell_dofs(iT) + local_boundary_dofs(iT));

    Eigen::MatrixXd RTN_mass_mat = Eigen::MatrixXd::Zero(RTN_basis.dimension(), RTN_basis.dimension());
    Eigen::MatrixXd RTN_RHS_mat = Eigen::MatrixXd::Zero(RTN_basis.dimension(), local_cell_dofs(iT) + local_boundary_dofs(iT));

    size_t offset = 0;

    if (m_edge_deg >= 1)
    {
        CellBasisType cell_test_space = construct_cell_basis(iT, m_edge_deg - 1);
        for (auto &func : laplace_cell_enrichment)
        {
        // auto the_bubble = bubble(m_mesh->cell(iT), 0, 1.0);
        // auto bubble_enrich = Functional::curl(the_bubble);
            cell_test_space.add_basis_function(func);
            // cell_test_space.add_basis_function(bubble_enrich);
        }

        if (laplace_cell_enrichment.size() != 0)
        {
            // Eigen::MatrixXd B_mat = Eigen::MatrixXd::Identity(cell_test_space.dimension(), cell_test_space.dimension());
            // B_mat.topLeftCorner(cell_test_space.dimension() - 1, cell_test_space.dimension() - 1) = cell_test_space.matrix();

            Eigen::MatrixXd mass = m_quad_handle.l2_product(cell_test_space.ancestor(), cell_test_space.ancestor(), m_mesh->cell(iT));
            Eigen::MatrixXd B_mat = orthonormalize(mass);

            cell_test_space.reset_matrix(B_mat);
        }

        RTN_mass_mat.topLeftCorner(cell_test_space.dimension(), RTN_basis.dimension()) = m_quad_handle.l2_product(cell_test_space, RTN_basis, m_mesh->cell(iT), false);
        RTN_RHS_mat.topLeftCorner(cell_test_space.dimension(), local_cell_dofs(iT)) = m_quad_handle.l2_product(cell_test_space, *m_cell_basis[iT].get(), m_mesh->cell(iT), false);

        offset = cell_test_space.dimension();
    }
    else
    {
        CellBasisType cell_test_space;
        for (auto &func : laplace_cell_enrichment)
        {
            cell_test_space.add_basis_function(func);
        }
        if (laplace_cell_enrichment.size() != 0)
        {
            cell_test_space.reset_matrix(Eigen::MatrixXd::Identity(1, 1));
        }

        // if(laplace_cell_enrichment.size() != 0 )
        // {
        //     Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(cell_test_space.ancestor(), cell_test_space.ancestor(), m_mesh->cell(iT), true);
        //     Eigen::MatrixXd B_mat = orthonormalize(basis_integrals);
        //     cell_test_space.reset_matrix(B_mat);
        // }

        RTN_mass_mat.topLeftCorner(cell_test_space.dimension(), RTN_basis.dimension()) = m_quad_handle.l2_product(cell_test_space, RTN_basis, m_mesh->cell(iT), false);
        RTN_RHS_mat.topLeftCorner(cell_test_space.dimension(), local_cell_dofs(iT)) = m_quad_handle.l2_product(cell_test_space, *m_cell_basis[iT].get(), m_mesh->cell(iT), false);

        offset = cell_test_space.dimension();
    }

    size_t RHS_col_offset = local_cell_dofs(iT);
    for (size_t iTE = 0; iTE < m_mesh->cell(iT)->n_edges(); ++iTE)
    {
        auto *edge = m_mesh->cell(iT)->edge(iTE);
        const size_t iE = edge->global_index();
        PolyMesh2D::Functional::ScalarFamily1D edge_test_space = construct_scalar_edge_basis(iE, m_edge_deg);

        double wTF = (double)m_mesh->cell(iT)->edge_orientation(iTE);

        {
            Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(edge_test_space.dimension(), RTN_basis.dimension());
            for (size_t j = 0; j < RTN_basis.dimension(); ++j)
            {
                auto rtn_j_dot_n(Functional::dot_n(Functional::trace(RTN_basis.ancestor().function(j), edge->parameterisation()), edge->parameterisation()));
                for (size_t i = 0; i < edge_test_space.dimension(); ++i)
                {
                    tmp(i, j) = m_quad_handle.integrate(edge_test_space.ancestor().function(i) * rtn_j_dot_n, edge);
                }
            }
            RTN_mass_mat.block(offset, 0, edge_test_space.dimension(), RTN_basis.dimension()) = wTF * edge_test_space.matrix() * tmp * RTN_basis.matrix().transpose();
        }

        {
            Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(edge_test_space.dimension(), local_edge_dofs(iE));
            for (size_t j = 0; j < local_edge_dofs(iE); ++j)
            {
                auto vF_j_dot_n(Functional::dot_n(m_edge_basis[iE]->ancestor().function(j), edge->parameterisation()));
                for (size_t i = 0; i < edge_test_space.dimension(); ++i)
                {
                    tmp(i, j) = m_quad_handle.integrate(edge_test_space.ancestor().function(i) * vF_j_dot_n, edge);
                }
            }
            RTN_RHS_mat.block(offset, RHS_col_offset, edge_test_space.dimension(), local_edge_dofs(iE)) = wTF * edge_test_space.matrix() * tmp * m_edge_basis[iE]->matrix().transpose();
        }

        // RTN_mass_mat.block(offset, 0, edge_test_space.dimension(), RTN_basis.dimension()) = m_quad_handle.l2_product(edge_test_space, RTN_basis, edge);
        // RTN_RHS_mat.block(offset, RHS_col_offset, edge_test_space.dimension(), local_edge_dofs(iE)) = m_quad_handle.l2_product(edge_test_space, *m_edge_basis[iE].get(), edge);

        offset += edge_test_space.dimension();
        RHS_col_offset += local_edge_dofs(iE);
    }

    // std::cout << RTN_mass_mat.norm() << "\n";

    // std::cout << "\n" << RTN_mass_mat << "\n";

    assert(offset == RTN_basis.dimension());

    Eigen::MatrixXd RTN_mat = RTN_mass_mat.inverse() * (RTN_RHS_mat);

    return RTN_mat.transpose() * m_quad_handle.l2_product(source, RTN_basis, m_mesh->cell(iT));

    // return Eigen::VectorXd::Zero(local_cell_dofs(iT) + local_boundary_dofs(iT));
}

Eigen::VectorXd StokesCore::l2_cell_projection(const size_t iT, const Functional::VectorFunction2D &func) const
{
    if (m_ortho)
    {
        return m_quad_handle.l2_product(func, *m_cell_basis[iT].get(), m_mesh->cell(iT));
    }

    Eigen::MatrixXd mass_mat = m_quad_handle.l2_product(*m_cell_basis[iT].get(), *m_cell_basis[iT].get(), m_mesh->cell(iT), true);
    return mass_mat.inverse() * m_quad_handle.l2_product(func, *m_cell_basis[iT].get(), m_mesh->cell(iT));
}

Eigen::VectorXd StokesCore::l2_edge_projection(const size_t iE, const Functional::VectorFunction1D &func) const
{
    if (m_ortho)
    {
        return m_quad_handle.l2_product(func, *m_edge_basis[iE].get(), m_mesh->edge(iE));
    }

    Eigen::MatrixXd mass_mat(m_quad_handle.l2_product(*m_edge_basis[iE].get(), *m_edge_basis[iE].get(), m_mesh->edge(iE), true));
    return mass_mat.inverse() * m_quad_handle.l2_product(func, *m_edge_basis[iE].get(), m_mesh->edge(iE));
}

Eigen::VectorXd StokesCore::l2_highorder_projection(const size_t iT, const Functional::VectorFunction2D &func) const
{
    if (m_ortho)
    {
        return m_quad_handle.l2_product(func, *m_highorder_basis[iT].get(), m_mesh->cell(iT));
    }

    Eigen::MatrixXd mass_mat(m_quad_handle.l2_product(*m_highorder_basis[iT].get(), *m_highorder_basis[iT].get(), m_mesh->cell(iT), true));
    return mass_mat.inverse() * m_quad_handle.l2_product(func, *m_highorder_basis[iT].get(), m_mesh->cell(iT));
}

Eigen::VectorXd StokesCore::l2_pressure_projection(const size_t iT, const Functional::ScalarFunction2D &func) const
{
    // if (m_ortho)
    // {
    //     return m_quad_handle.l2_product(func, *m_pressure_basis[iT].get(), m_mesh->cell(iT));
    // }

    Eigen::MatrixXd mass_mat(m_quad_handle.l2_product(*m_pressure_basis[iT].get(), *m_pressure_basis[iT].get(), m_mesh->cell(iT), true));
    return mass_mat.inverse() * m_quad_handle.l2_product(func, *m_pressure_basis[iT].get(), m_mesh->cell(iT));
}

Eigen::VectorXd StokesCore::elliptic_projection(const size_t iT, const Functional::VectorFunction2D &func) const
{
    Eigen::MatrixXd ST(m_quad_handle.h1_product(*m_highorder_basis[iT].get(), *m_highorder_basis[iT].get(), m_mesh->cell(iT), true));

    Eigen::MatrixXd LT = Eigen::MatrixXd::Zero(local_highorder_dofs(iT), local_highorder_dofs(iT));
    Eigen::VectorXd RHS = Eigen::VectorXd::Zero(local_highorder_dofs(iT));

    double scalT = ST.norm();

    for (size_t i = 0; i < local_highorder_dofs(iT); ++i)
    {
        RHS(i) = m_quad_handle.integrate(Functional::gradient_product(func, m_highorder_basis[iT]->ancestor().function(i)), m_mesh->cell(iT));
        RHS(i) += scalT * Math::scalar_product(m_quad_handle.integrate(func, m_mesh->cell(iT)), m_quad_handle.integrate(m_highorder_basis[iT]->ancestor().function(i), m_mesh->cell(iT)));
        for (size_t j = i; j < local_highorder_dofs(iT); ++j)
        {
            LT(j, i) = LT(i, j) = Math::scalar_product(m_quad_handle.integrate(m_highorder_basis[iT]->ancestor().function(i), m_mesh->cell(iT)), m_quad_handle.integrate(m_highorder_basis[iT]->ancestor().function(j), m_mesh->cell(iT)));
        }
    }

    RHS = m_highorder_basis[iT]->matrix() * RHS;
    LT = m_highorder_basis[iT]->matrix() * LT * m_highorder_basis[iT]->matrix().transpose();

    return (ST + scalT * LT).ldlt().solve(RHS);
}

void StokesCore::enrich_highorder_basis(const size_t iT, const Functional::VectorFunction2D &func)
{
    assert(iT < m_highorder_basis.size());

    // assuming orthonormalization
    m_highorder_basis[iT]->add_basis_function(func);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(m_highorder_basis[iT]->ancestor(), m_highorder_basis[iT]->ancestor(), m_mesh->cell(iT), true);

    Eigen::VectorXcd eigs = basis_integrals.eigenvalues();
    double max = std::abs(eigs(0));
    double min = std::abs(eigs(0));
    for (int i = 1; i < eigs.size(); ++i)
    {
        max = std::max(max, std::abs(eigs(i)));
        min = std::min(min, std::abs(eigs(i)));
    }

    if ((max / min) > 1E16)
    {
        std::cerr << "Highorder basis #" << iT << " is linearly dependent.\n";
        exit(1);
    }

    Eigen::MatrixXd B = orthonormalize(basis_integrals);
    m_highorder_basis[iT]->reset_matrix(B);
}

void StokesCore::enrich_cell_basis(const size_t iT, const Functional::VectorFunction2D &func)
{
    assert(iT < m_cell_basis.size());

    // assuming orthonormalization
    m_cell_basis[iT]->add_basis_function(func);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(m_cell_basis[iT]->ancestor(), m_cell_basis[iT]->ancestor(), m_mesh->cell(iT), true);

    Eigen::VectorXcd eigs = basis_integrals.eigenvalues();
    double max = std::abs(eigs(0));
    double min = std::abs(eigs(0));
    for (int i = 1; i < eigs.size(); ++i)
    {
        max = std::max(max, std::abs(eigs(i)));
        min = std::min(min, std::abs(eigs(i)));
    }

    if ((max / min) > 1E16)
    {
        std::cerr << "Cell basis #" << iT << " is linearly dependent.\n";
        exit(1);
    }

    Eigen::MatrixXd B = orthonormalize(basis_integrals);
    m_cell_basis[iT]->reset_matrix(B);
}

void StokesCore::enrich_edge_basis(const size_t iE, const Functional::VectorFunction1D &func)
{
    assert(iE < m_edge_basis.size());

    // assuming orthonormalization
    m_edge_basis[iE]->add_basis_function(func);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(m_edge_basis[iE]->ancestor(), m_edge_basis[iE]->ancestor(), m_mesh->edge(iE), true);

    Eigen::VectorXcd eigs = basis_integrals.eigenvalues();
    double max = std::abs(eigs(0));
    double min = std::abs(eigs(0));
    for (int i = 1; i < eigs.size(); ++i)
    {
        max = std::max(max, std::abs(eigs(i)));
        min = std::min(min, std::abs(eigs(i)));
    }

    if ((max / min) > 1E16)
    {
        std::cerr << "Edge basis #" << iE << " is linearly dependent.\n";
        exit(1);
    }

    Eigen::MatrixXd B = orthonormalize(basis_integrals);
    m_edge_basis[iE]->reset_matrix(B);
}

void StokesCore::enrich_pressure_basis(const size_t iT, const Functional::ScalarFunction2D &func)
{
    assert(iT < m_pressure_basis.size());

    // assuming orthonormalization
    m_pressure_basis[iT]->add_basis_function(func);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(m_pressure_basis[iT]->ancestor(), m_pressure_basis[iT]->ancestor(), m_mesh->cell(iT), true);

    Eigen::VectorXcd eigs = basis_integrals.eigenvalues();
    double max = std::abs(eigs(0));
    double min = std::abs(eigs(0));
    for (int i = 1; i < eigs.size(); ++i)
    {
        max = std::max(max, std::abs(eigs(i)));
        min = std::min(min, std::abs(eigs(i)));
    }

    if ((max / min) > 1E16)
    {
        std::cerr << "Pressure basis #" << iT << " is linearly dependent.\n";
        std::cout << "Condition number of gram matrix = " << max / min << "\n";
        exit(1);
    }

    Eigen::MatrixXd B = orthonormalize(basis_integrals);
    m_pressure_basis[iT]->reset_matrix(B);

    // Eigen::MatrixXd basis_integrals_new = m_quad_handle.l2_product(*m_pressure_basis[iT].get(), *m_pressure_basis[iT].get(), m_mesh->cell(iT), true);

    // std::cout << "\n\n" << basis_integrals_new << "\n\n";

    // exit(1);
}

Eigen::VectorXd StokesCore::velocity_restr(const Eigen::VectorXd &UVec, const size_t iT) const
{
    assert((size_t)UVec.size() == total_cell_dofs() + total_edge_dofs() + total_pressure_dofs());

    Eigen::VectorXd UVec_T = Eigen::VectorXd::Zero(local_cell_dofs(iT) + local_boundary_dofs(iT));

    UVec_T.head(local_cell_dofs(iT)) = UVec.segment(global_offset_T(iT), local_cell_dofs(iT));

    for (size_t iTE = 0; iTE < m_mesh->cell(iT)->n_edges(); iTE++)
    {
        const size_t iE = m_mesh->cell(iT)->edge(iTE)->global_index();
        UVec_T.segment(local_cell_dofs(iT) + local_offset_E(iT, iTE), local_edge_dofs(iE)) = UVec.segment(total_cell_dofs() + global_offset_E(iE), local_edge_dofs(iE));
    }

    return UVec_T;
}

Eigen::VectorXd StokesCore::pressure_restr(const Eigen::VectorXd &UVec, const size_t iT) const
{
    assert((size_t)UVec.size() == total_cell_dofs() + total_edge_dofs() + total_pressure_dofs());
    return UVec.segment(total_cell_dofs() + total_edge_dofs() + global_pressure_offset_T(iT), local_pressure_dofs(iT));
}

Eigen::VectorXd StokesCore::interpolate(const Functional::Function<2, 2> &u_func, const Functional::Function<2, 1> &p_func) const
{
    Eigen::VectorXd Ihk(Eigen::VectorXd::Zero(total_cell_dofs() + total_edge_dofs() + total_pressure_dofs()));

    std::function<void(size_t, size_t)> interpolate_cell_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            Ihk.segment(global_offset_T(iT), local_cell_dofs(iT)) = this->l2_cell_projection(iT, u_func);
            Ihk.segment(total_cell_dofs() + total_edge_dofs() + global_pressure_offset_T(iT), local_pressure_dofs(iT)) = this->l2_pressure_projection(iT, p_func);
            // Eigen::MatrixXd MTT = m_quad_handle.l2_product(*m_cell_basis[iT].get(), *m_cell_basis[iT].get(), m_mesh->cell(iT), true);
            // Eigen::VectorXd u_integrated_against_basis_T = m_quad_handle.l2_product(u_func, *m_cell_basis[iT].get(), m_mesh->cell(iT));
            // Eigen::VectorXd p_integrated_against_pressure_basis = m_quad_handle.l2_product(p_func, *m_pressure_basis[iT].get(), m_mesh->cell(iT));

            // Eigen::MatrixXd pressure_mass_mat = m_quad_handle.l2_product(*m_pressure_basis[iT].get(), *m_pressure_basis[iT].get(), m_mesh->cell(iT), true);

            // Ihk.segment(global_offset_T(iT), local_cell_dofs(iT)) = MTT.inverse() * u_integrated_against_basis_T;
            // Ihk.segment(total_cell_dofs() + total_edge_dofs() + global_pressure_offset_T(iT), local_pressure_dofs(iT)) = pressure_mass_mat.inverse() * p_integrated_against_pressure_basis;
        }
    };
    parallel_for(m_mesh->n_cells(), interpolate_cell_terms, m_use_threads);

    std::function<void(size_t, size_t)> interpolate_edge_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iF = start; iF < end; ++iF)
        {
            Ihk.segment(total_cell_dofs() + global_offset_E(iF), local_edge_dofs(iF)) = this->l2_edge_projection(iF, Functional::trace(u_func, m_mesh->edge(iF)->parameterisation()));
            // Eigen::MatrixXd MFF = m_quad_handle.l2_product(*m_edge_basis[iF].get(), *m_edge_basis[iF].get(), m_mesh->edge(iF), true);
            // Eigen::VectorXd u_integrated_against_basis_F = m_quad_handle.l2_product(Functional::trace(u_func, m_mesh->edge(iF)->parameterisation()), *m_edge_basis[iF].get(), m_mesh->edge(iF));

            // Ihk.segment(total_cell_dofs() + global_offset_E(iF), local_edge_dofs(iF)) = MFF.inverse() * u_integrated_against_basis_F;
        }
    };
    parallel_for(m_mesh->n_edges(), interpolate_edge_terms, m_use_threads);

    return Ihk;
}

const Quadrature::QuadHandler<CurvedMesh::Mesh> &StokesCore::get_quad_handle() const
{
    return m_quad_handle;
}

StokesCore::MeshType *StokesCore::get_mesh() const
{
    return m_mesh;
}
size_t StokesCore::cell_degree() const
{
    return m_cell_deg;
}
size_t StokesCore::edge_degree() const
{
    return m_edge_deg;
}

StokesCore::CellBasisType *StokesCore::highorder_basis(const size_t iT) const
{
    assert(iT < m_highorder_basis.size());
    return m_highorder_basis[iT].get();
}
StokesCore::CellBasisType *StokesCore::cell_basis(const size_t iT) const
{
    assert(iT < m_cell_basis.size());
    return m_cell_basis[iT].get();
}
StokesCore::EdgeBasisType *StokesCore::edge_basis(const size_t iE) const
{
    assert(iE < m_edge_basis.size());
    return m_edge_basis[iE].get();
}
StokesCore::ScalarBasisType *StokesCore::pressure_basis(const size_t iT) const
{
    assert(iT < m_pressure_basis.size());
    return m_pressure_basis[iT].get();
}

size_t StokesCore::local_cell_dofs(const size_t iT) const
{
    return m_cell_basis[iT]->dimension();
}
size_t StokesCore::local_highorder_dofs(const size_t iT) const
{
    return m_highorder_basis[iT]->dimension();
}
size_t StokesCore::local_edge_dofs(const size_t iE) const
{
    return m_edge_basis[iE]->dimension();
}
size_t StokesCore::local_pressure_dofs(const size_t iT) const
{
    return m_pressure_basis[iT]->dimension();
}

size_t StokesCore::local_boundary_dofs(const size_t iT) const
{
    size_t dofs = 0;
    for (size_t iTF = 0; iTF < m_mesh->cell(iT)->n_edges(); ++iTF)
    {
        dofs += this->local_edge_dofs(m_mesh->cell(iT)->edge(iTF)->global_index());
    }
    return dofs;
}

size_t StokesCore::global_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    size_t offset = 0;
    for (size_t i = 0; i < iT; ++i)
    {
        offset += this->local_cell_dofs(i);
    }
    return offset;
}
size_t StokesCore::global_offset_E(const size_t iE) const
{
    assert(iE < m_mesh->n_edges());
    size_t offset = 0;
    for (size_t i = 0; i < iE; ++i)
    {
        offset += this->local_edge_dofs(i);
    }
    // offset += this->total_cell_dofs();
    return offset;
}
size_t StokesCore::local_offset_E(const size_t iT, const size_t iTE) const
{
    assert(iTE < m_mesh->cell(iT)->n_edges());
    size_t offset = 0;
    for (size_t i = 0; i < iTE; ++i)
    {
        size_t iE = m_mesh->cell(iT)->edge(i)->global_index();
        offset += this->local_edge_dofs(iE);
    }
    return offset;
}
size_t StokesCore::local_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    return 0;
}

size_t StokesCore::global_pressure_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    size_t offset = 0;
    for (size_t i = 0; i < iT; ++i)
    {
        offset += this->local_pressure_dofs(i);
    }
    return offset;
}

size_t StokesCore::local_pressure_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    return 0;
}

size_t StokesCore::total_cell_dofs() const
{
    size_t dofs = 0;
    for (size_t iT = 0; iT < m_mesh->n_cells(); ++iT)
    {
        dofs += this->local_cell_dofs(iT);
    }
    return dofs;
}
size_t StokesCore::total_edge_dofs() const
{
    size_t dofs = 0;
    for (size_t iE = 0; iE < m_mesh->n_edges(); ++iE)
    {
        dofs += this->local_edge_dofs(iE);
    }
    return dofs;
}
size_t StokesCore::total_pressure_dofs() const
{
    size_t dofs = 0;
    for (size_t iT = 0; iT < m_mesh->n_cells(); ++iT)
    {
        dofs += this->local_pressure_dofs(iT);
    }
    return dofs;
}
