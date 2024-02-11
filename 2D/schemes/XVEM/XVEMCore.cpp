#include <XVEMCore.h>
#include <parallel_for.hpp>
#include "common.hpp"
#include "vtu_writer.hpp"
#include <iomanip>
#include <boost/timer/timer.hpp>
#include <algorithm>

using namespace PolyMesh2D;



void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}

Eigen::MatrixXd orthonormalize(const Eigen::MatrixXd &gram_mat)
{
    Eigen::LLT<Eigen::MatrixXd> llt(gram_mat);
    Eigen::MatrixXd L = llt.matrixL();
    return L.inverse();
}

XVEMCore::XVEMCore(
    MeshType *mesh_ptr, ///< A pointer to the loaded mesh
    size_t l,    ///< The degree of the cell polynomials
    size_t k,    ///< The degree of the edge polynomials
    bool use_threads,   ///< Optional argument to indicate if threads should be used
    bool ortho)

    : m_mesh(mesh_ptr),
      m_l(l),
      m_k(k),
    //   m_quad_handle(m_mesh, 2 * m_k + 3, 2 * m_k + 3),
      m_quad_handle(m_mesh, 25, 25),
      m_use_threads(use_threads),
      m_ortho(ortho),
      m_high_cell_basis(mesh_ptr->n_cells()),
      m_low_cell_basis(mesh_ptr->n_cells()),
      m_high_edge_basis(mesh_ptr->n_edges()),
      m_low_edge_basis(mesh_ptr->n_edges())
{
    boost::timer::cpu_timer construction_timer;
    construction_timer.start();

    std::cout << "\n[XVEMCore] Constructing bases\n";

    // Create extended bases
    std::function<void(size_t, size_t)> construct_all_cell_bases = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            m_high_cell_basis[iT].reset(new CellBasisType(construct_cell_basis(iT, m_k)));
            m_low_cell_basis[iT].reset(new CellBasisType(construct_cell_basis(iT, m_l)));
        }
    };
    parallel_for(m_mesh->n_cells(), construct_all_cell_bases, m_use_threads);

    // Create cell bases
    std::function<void(size_t, size_t)> construct_all_edge_bases = [&](size_t start, size_t end) -> void
    {
        for (size_t iE = start; iE < end; iE++)
        {
            m_high_edge_basis[iE].reset(new EdgeBasisType(construct_edge_basis(iE, m_k)));
            m_low_edge_basis[iE] = std::make_unique<EdgeBasisType>(*m_high_edge_basis[iE]); // copy of high edge basis

            Eigen::MatrixXd mat(m_low_edge_basis[iE]->matrix());
            // careful.... must remove rows from bottom up. 
            // i.e. removeRow(mat, m_k - 1); removeRow(mat, m_k); will actually remove rows k-1, k+1
            removeRow(mat, m_k); // remove row associated with basis P^k
            removeRow(mat, m_k - 1); // remove row associated with basis P^{k-1}
            m_low_edge_basis[iE]->reset_matrix(mat);
        }
    };
    parallel_for(m_mesh->n_edges(), construct_all_edge_bases, m_use_threads);

    construction_timer.stop();
    std::cout << "     Construction time = " << construction_timer.elapsed().wall * 1E-9 << "s\n";
}

XVEMCore::CellBasisType XVEMCore::construct_cell_basis(const size_t iT, const size_t degree) const
{
    // Transform to cell centred coodinates for better conditioning
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

    Eigen::MatrixXd B = (true ? orthonormalize(m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->cell(iT), true)) : Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension()));

    return CellBasisType(mono_basis, B);
}

XVEMCore::EdgeBasisType XVEMCore::construct_edge_basis(const size_t iE, const int degree) const
{
    // Basis of monomials

    if(degree < 0)
    {
        assert(degree == -1);
        return EdgeBasisType(); // return empty basis
    }

    PolyMesh2D::Functional::ScalarBasis1D mono_basis = PolyMesh2D::Functional::MonomialScalarBasis(m_mesh->edge(iE)->parameterisation(), degree);

    Eigen::MatrixXd B = (true ? orthonormalize(m_quad_handle.l2_product(mono_basis, mono_basis, m_mesh->edge(iE), true)) : Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension()));

    return EdgeBasisType(mono_basis, B);
}

Eigen::VectorXd XVEMCore::l2_low_cell_projection(const size_t iT, const Functional::ScalarFunction2D &func) const
{
    if (m_ortho)
    {
        return m_quad_handle.l2_product(func, *m_low_cell_basis[iT].get(), m_mesh->cell(iT));
    }

    Eigen::MatrixXd mass_mat = m_quad_handle.l2_product(*m_low_cell_basis[iT].get(), *m_low_cell_basis[iT].get(), m_mesh->cell(iT), true);
    return mass_mat.inverse() * m_quad_handle.l2_product(func, *m_low_cell_basis[iT].get(), m_mesh->cell(iT));
}

Eigen::VectorXd XVEMCore::l2_low_edge_projection(const size_t iE, const Functional::ScalarFunction1D &func) const
{
    if (m_ortho)
    {
        return m_quad_handle.l2_product(func, *m_low_edge_basis[iE].get(), m_mesh->edge(iE));
    }

    Eigen::MatrixXd mass_mat(m_quad_handle.l2_product(*m_low_edge_basis[iE].get(), *m_low_edge_basis[iE].get(), m_mesh->edge(iE), true));
    return mass_mat.inverse() * m_quad_handle.l2_product(func, *m_low_edge_basis[iE].get(), m_mesh->edge(iE));
}

Eigen::VectorXd XVEMCore::elliptic_projection(const size_t iT, const Functional::ScalarFunction2D &func) const
{
    Eigen::MatrixXd ST(m_quad_handle.h1_product(*m_high_cell_basis[iT].get(), *m_high_cell_basis[iT].get(), m_mesh->cell(iT), true));

    Eigen::MatrixXd LT = Eigen::MatrixXd::Zero(local_high_cell_dofs(iT), local_high_cell_dofs(iT));
    Eigen::VectorXd RHS = Eigen::VectorXd::Zero(local_high_cell_dofs(iT));

    double scalT = ST.norm();

    for (size_t i = 0; i < local_high_cell_dofs(iT); ++i)
    {
        RHS(i) = m_quad_handle.integrate(Functional::gradient_product(func, m_high_cell_basis[iT]->ancestor().function(i)), m_mesh->cell(iT));
        RHS(i) += scalT * Math::scalar_product(m_quad_handle.integrate(func, m_mesh->cell(iT)), m_quad_handle.integrate(m_high_cell_basis[iT]->ancestor().function(i), m_mesh->cell(iT)));
        for (size_t j = i; j < local_high_cell_dofs(iT); ++j)
        {
            LT(j, i) = LT(i, j) = Math::scalar_product(m_quad_handle.integrate(m_high_cell_basis[iT]->ancestor().function(i), m_mesh->cell(iT)), m_quad_handle.integrate(m_high_cell_basis[iT]->ancestor().function(j), m_mesh->cell(iT)));
        }
    }

    RHS = m_high_cell_basis[iT]->matrix() * RHS;
    LT = m_high_cell_basis[iT]->matrix() * LT * m_high_cell_basis[iT]->matrix().transpose();

    return (ST + scalT * LT).ldlt().solve(RHS);
}

void XVEMCore::enrich_cell_basis(const size_t iT, const Functional::ScalarFunction2D &func)
{
    assert(iT < m_high_cell_basis.size());

    // assume zero Laplacian so low_cell_basis not enriched.

    // assuming orthonormalization
    m_high_cell_basis[iT]->add_basis_function(func);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(m_high_cell_basis[iT]->ancestor(), m_high_cell_basis[iT]->ancestor(), m_mesh->cell(iT), true);

    // Eigen::VectorXcd eigs = basis_integrals.eigenvalues();
    // double max = std::abs(eigs(0));
    // double min = std::abs(eigs(0));
    // for (int i = 1; i < eigs.size(); ++i)
    // {
    //     max = std::max(max, std::abs(eigs(i)));
    //     min = std::min(min, std::abs(eigs(i)));
    // }

    // if ((max / min) > 1E16)
    // {
    //     std::cerr << "Cell basis #" << iT << " is linearly dependent.\n";
    //     // exit(1);

    //     m_high_cell_basis[iT]->remove_basis_function(m_high_cell_basis[iT]->dimension() - 1);
    //     basis_integrals = m_quad_handle.l2_product(m_high_cell_basis[iT]->ancestor(), m_high_cell_basis[iT]->ancestor(), m_mesh->cell(iT), true);
    // }

    Eigen::MatrixXd B = orthonormalize(basis_integrals);
    m_high_cell_basis[iT]->reset_matrix(B);

    // assume zero Laplacian so low_cell_basis not enriched.

    // assuming orthonormalization
    // m_low_cell_basis[iT]->add_basis_function(func);
    // basis_integrals = m_quad_handle.l2_product(m_low_cell_basis[iT]->ancestor(), m_low_cell_basis[iT]->ancestor(), m_mesh->cell(iT), true);

    // B = orthonormalize(basis_integrals);
    // m_low_cell_basis[iT]->reset_matrix(B);
}

void XVEMCore::enrich_edge_basis(const size_t iE, const Functional::ScalarFunction1D &func)
{
    assert(iE < m_high_edge_basis.size());

    m_high_edge_basis[iE]->add_basis_function(func);
    Eigen::MatrixXd basis_integrals = m_quad_handle.l2_product(m_high_edge_basis[iE]->ancestor(), m_high_edge_basis[iE]->ancestor(), m_mesh->edge(iE), true);

    Eigen::VectorXcd eigs = basis_integrals.eigenvalues();
    double max = std::abs(eigs(0));
    double min = std::abs(eigs(0));
    for (int i = 1; i < eigs.size(); ++i)
    {
        max = std::max(max, std::abs(eigs(i)));
        min = std::min(min, std::abs(eigs(i)));
    }

    if (!((max / min) < 1E16))
    {
        // std::cout << " [XVEMCore::enrich_edge_basis] Edge basis #" << iE << " is linearly dependent (cond = " << max / min <<"). Removing enrichment\n";

        m_high_edge_basis[iE]->remove_basis_function(m_high_edge_basis[iE]->ancestor().dimension() - 1);
        basis_integrals = m_quad_handle.l2_product(m_high_edge_basis[iE]->ancestor(), m_high_edge_basis[iE]->ancestor(), m_mesh->edge(iE), true);

        Eigen::MatrixXd B = orthonormalize(basis_integrals);
        m_high_edge_basis[iE]->reset_matrix(B);
        return;
    }

    Eigen::MatrixXd B = orthonormalize(basis_integrals);
    m_high_edge_basis[iE]->reset_matrix(B);

    m_low_edge_basis[iE]->add_basis_function(func);
    removeRow(B, m_k); // remove row associated with basis P^k
    removeRow(B, m_k - 1); // remove row associated with basis P^{k-1}

    m_low_edge_basis[iE]->reset_matrix(B);
    
    // m_low_edge_basis[iE]->add_basis_function(func);
    // basis_integrals = m_quad_handle.l2_product(m_low_edge_basis[iE]->ancestor(), m_low_edge_basis[iE]->ancestor(), m_mesh->edge(iE), true);

    // B = orthonormalize(basis_integrals);
    // m_low_edge_basis[iE]->reset_matrix(B);
}

Eigen::VectorXd XVEMCore::restr(const Eigen::VectorXd &UVec, const size_t iT) const
{
    assert((size_t)UVec.size() == total_cell_dofs() + total_edge_dofs() + m_mesh->n_vertices());

    Eigen::VectorXd UVec_T = Eigen::VectorXd::Zero(local_dofs(iT));

    UVec_T.head(local_low_cell_dofs(iT)) = UVec.segment(global_offset_T(iT), local_low_cell_dofs(iT));

    for (size_t iTE = 0; iTE < m_mesh->cell(iT)->n_edges(); iTE++)
    {
        const size_t iE = m_mesh->cell(iT)->edge(iTE)->global_index();
        UVec_T.segment(local_offset_E(iT, iTE), local_low_edge_dofs(iE)) = UVec.segment(global_offset_E(iE), local_low_edge_dofs(iE));
    }

    for (size_t iTV = 0; iTV < m_mesh->cell(iT)->n_vertices(); iTV++)
    {
        const size_t iV = m_mesh->cell(iT)->vertex(iTV)->global_index();
        UVec_T(local_offset_V(iT, iTV)) = UVec(global_offset_V(iV));
    }

    return UVec_T;
}

Eigen::VectorXd XVEMCore::interpolate(const Functional::Function<2, 1> &func) const
{
    Eigen::VectorXd Ihk(Eigen::VectorXd::Zero(total_cell_dofs() + total_edge_dofs() + m_mesh->n_vertices()));

    std::function<void(size_t, size_t)> interpolate_cell_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            Ihk.segment(global_offset_T(iT), local_low_cell_dofs(iT)) = this->l2_low_cell_projection(iT, func);
        }
    };
    parallel_for(m_mesh->n_cells(), interpolate_cell_terms, m_use_threads);

    std::function<void(size_t, size_t)> interpolate_edge_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iE = start; iE < end; ++iE)
        {
            Ihk.segment(global_offset_E(iE), local_low_edge_dofs(iE)) = this->l2_low_edge_projection(iE, Functional::trace(func, m_mesh->edge(iE)->parameterisation()));
        }
    };
    parallel_for(m_mesh->n_edges(), interpolate_edge_terms, m_use_threads);

    std::function<void(size_t, size_t)> interpolate_vertex_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iV = start; iV < end; ++iV)
        {
            Ihk(global_offset_V(iV)) = func.value(m_mesh->vertex(iV)->coords());
        }
    };
    parallel_for(m_mesh->n_vertices(), interpolate_vertex_terms, m_use_threads);

    return Ihk;
}

Eigen::MatrixXd XVEMCore::continuous_reconstruction(const size_t iT, const size_t iTE) const
{
    const size_t iE = m_mesh->cell(iT)->edge(iTE)->global_index();
    const size_t n_local_low_edge_dofs = local_low_edge_dofs(iE);    
    const size_t n_local_high_edge_dofs = local_high_edge_dofs(iE);       

    assert(n_local_low_edge_dofs == n_local_high_edge_dofs - 2);

    Eigen::MatrixXd LHS_mat = Eigen::MatrixXd::Zero(n_local_high_edge_dofs, n_local_high_edge_dofs);
    Eigen::MatrixXd RHS_mat = Eigen::MatrixXd::Identity(n_local_high_edge_dofs, n_local_high_edge_dofs);

    LHS_mat.topLeftCorner(n_local_low_edge_dofs, n_local_high_edge_dofs) = m_quad_handle.l2_product(*m_low_edge_basis[iE].get(), *m_high_edge_basis[iE].get(), m_mesh->edge(iE), false);
    RHS_mat.topLeftCorner(n_local_low_edge_dofs, n_local_low_edge_dofs)  = m_quad_handle.l2_product(*m_low_edge_basis[iE].get(),  *m_low_edge_basis[iE].get(), m_mesh->edge(iE), true);

    PolyMesh2D::Functional::Curve edge_param = m_mesh->edge(iE)->parameterisation();

    double t0 = edge_param.tmin;
    double t1 = edge_param.tmax;

    if((edge_param.value(t0) - m_mesh->edge(iE)->vertex(0)->coords()).norm() < 1E-14)
    {
        assert((edge_param.value(t1) - m_mesh->edge(iE)->vertex(1)->coords()).norm() < 1E-14);
    }
    else
    {
        std::swap(t0, t1);

        assert((edge_param.value(t0) - m_mesh->edge(iE)->vertex(0)->coords()).norm() < 1E-14);
        assert((edge_param.value(t1) - m_mesh->edge(iE)->vertex(1)->coords()).norm() < 1E-14);
    }

    for(size_t i = 0; i < n_local_high_edge_dofs; ++i)
    {
        LHS_mat(n_local_low_edge_dofs, i)     = m_high_edge_basis[iE]->value(i, t0);
        LHS_mat(n_local_low_edge_dofs + 1, i) = m_high_edge_basis[iE]->value(i, t1);
    }

    return LHS_mat.inverse() * RHS_mat * this->edge_restrict_mat(iT, iTE); // TODO: use a better solver
}

Eigen::MatrixXd XVEMCore::edge_restrict_mat(const size_t iT, const size_t iTE) const
{
    const size_t iE = m_mesh->cell(iT)->edge(iTE)->global_index();
    const size_t n_local_low_edge_dofs = local_low_edge_dofs(iE);

    assert(local_low_edge_dofs(iE) == local_high_edge_dofs(iE) - 2);

    const size_t iV1 = m_mesh->edge(iE)->vertex(0)->global_index();
    const size_t iV2 = m_mesh->edge(iE)->vertex(1)->global_index();

    size_t iTV1 = 0;
    size_t iTV2 = 0;

    for(size_t iTV = 0; iTV < m_mesh->cell(iT)->n_vertices(); ++iTV)
    {
        if(m_mesh->cell(iT)->vertex(iTV) == m_mesh->vertex(iV1))
        {
            iTV1 = iTV;
        }
        else if(m_mesh->cell(iT)->vertex(iTV) == m_mesh->vertex(iV2))
        {
            iTV2 = iTV;
        }
    }
    assert((m_mesh->cell(iT)->vertex(iTV1) == m_mesh->vertex(iV1)) && (m_mesh->cell(iT)->vertex(iTV2) == m_mesh->vertex(iV2)));  

    const size_t local_offset = local_offset_E(iT, iTE);
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n_local_low_edge_dofs + 2, local_dofs(iT));
    for(size_t i = 0; i < n_local_low_edge_dofs; ++i)
    {
        mat(i, local_offset + i) = 1.0;
    }

    mat(n_local_low_edge_dofs, local_offset_V(iT, iTV1)) = 1.0;
    mat(n_local_low_edge_dofs + 1, local_offset_V(iT, iTV2)) = 1.0;

    return mat;
}

Eigen::MatrixXd XVEMCore::potential_reconstruction(size_t iT) const
{
    CurvedMesh::Cell * cell = m_mesh->cell(iT);

    const size_t n_local_low_cell_dofs = local_low_cell_dofs(iT);
    const size_t n_local_high_cell_dofs = local_high_cell_dofs(iT);
    const size_t n_local_dofs = local_dofs(iT);

    // Calculate stiffness matrix
    Eigen::MatrixXd ST(m_quad_handle.h1_product(*m_high_cell_basis[iT].get(), *m_high_cell_basis[iT].get(), cell, true));

    Eigen::MatrixXd M_LOWT_LOWT  = m_quad_handle.l2_product(*m_low_cell_basis[iT].get(),  *m_low_cell_basis[iT].get(), cell, true);
    Eigen::MatrixXd M_LOWT_HIGHT = m_quad_handle.l2_product(*m_low_cell_basis[iT].get(), *m_high_cell_basis[iT].get(), cell, false);

    Eigen::MatrixXd PT_RHS = Eigen::MatrixXd::Zero(n_local_high_cell_dofs, n_local_dofs);

    Eigen::MatrixXd L_HIGHT_LOWT  = M_LOWT_HIGHT.row(0).transpose() * M_LOWT_LOWT.row(0);
    Eigen::MatrixXd L_HIGHT_HIGHT = M_LOWT_HIGHT.row(0).transpose() * M_LOWT_HIGHT.row(0);

    double scalT = ST.norm() / L_HIGHT_HIGHT.norm();

    PT_RHS.topLeftCorner(n_local_high_cell_dofs, n_local_low_cell_dofs) = m_quad_handle.h1_product(*m_high_cell_basis[iT].get(), *m_low_cell_basis[iT].get(), cell, false) + scalT * L_HIGHT_LOWT;

    for (size_t iTE = 0; iTE < cell->n_edges(); iTE++)
    {
        CurvedMesh::Edge * edge = cell->edge(iTE);

        const size_t iE = edge->global_index();
        const size_t n_local_high_edge_dofs = local_high_edge_dofs(iE);      

        PolyMesh2D::Functional::Curve edge_param = edge->parameterisation();

        Eigen::MatrixXd gradw_nTF_vF = Eigen::MatrixXd::Zero(n_local_high_cell_dofs, n_local_high_edge_dofs);
        Eigen::MatrixXd gradw_nTF_vT = Eigen::MatrixXd::Zero(n_local_high_cell_dofs, n_local_low_cell_dofs);

        Eigen::Vector2d bias = -1E-15 * cell->edge_orientation(iTE) * edge_param.normal(0.5 * (edge_param.tmin + edge_param.tmax));
        
        for (size_t i = 0; i < n_local_high_cell_dofs; ++i)
        {
            Functional::ScalarFunction1D gradw_nTF;
            if(edge->is_boundary() && i == n_local_high_cell_dofs - 1 && std::abs(edge->center_mass()(1)) < 1E-14)
            {
                gradw_nTF = Functional::neumann_trace(m_high_cell_basis[iT]->ancestor().function(i), edge_param, bias);
            }
            else
            {
                gradw_nTF = Functional::neumann_trace(m_high_cell_basis[iT]->ancestor().function(i), edge_param);
            }
            for (size_t j = 0; j < n_local_high_edge_dofs; ++j)
            {
                gradw_nTF_vF(i, j) = m_quad_handle.integrate(Functional::scalar_product(gradw_nTF, m_high_edge_basis[iE]->ancestor().function(j)), edge);
            }
            for (size_t j = 0; j < n_local_low_cell_dofs; ++j)
            {
                auto vT_on_F(Functional::trace(m_low_cell_basis[iT]->ancestor().function(j), edge_param));
                gradw_nTF_vT(i, j) = m_quad_handle.integrate(Functional::scalar_product(gradw_nTF, vT_on_F), edge);
            }
        }

        PT_RHS.topLeftCorner(n_local_high_cell_dofs, n_local_low_cell_dofs) -= cell->edge_orientation(iTE) * m_high_cell_basis[iT]->matrix() * gradw_nTF_vT * m_low_cell_basis[iT]->matrix().transpose();

        PT_RHS += cell->edge_orientation(iTE) * m_high_cell_basis[iT]->matrix() * gradw_nTF_vF * m_high_edge_basis[iE]->matrix().transpose() * this->continuous_reconstruction(iT, iTE);
    }

    return (ST + scalT * L_HIGHT_HIGHT).ldlt().solve(PT_RHS);
}

const Quadrature::QuadHandler<CurvedMesh::Mesh> &XVEMCore::get_quad_handle() const
{
    return m_quad_handle;
}

XVEMCore::MeshType *XVEMCore::get_mesh() const
{
    return m_mesh;
}
size_t XVEMCore::get_l() const
{
    return m_l;
}
size_t XVEMCore::get_k() const
{
    return m_k;
}

XVEMCore::CellBasisType *XVEMCore::high_cell_basis(const size_t iT) const
{
    return m_high_cell_basis.at(iT).get();
}

XVEMCore::CellBasisType *XVEMCore::low_cell_basis(const size_t iT) const
{
    return m_low_cell_basis.at(iT).get();
}

XVEMCore::EdgeBasisType *XVEMCore::high_edge_basis(const size_t iE) const
{
    return m_high_edge_basis.at(iE).get();
}

XVEMCore::EdgeBasisType *XVEMCore::low_edge_basis(const size_t iE) const
{
    return m_low_edge_basis.at(iE).get();
}

size_t XVEMCore::local_high_cell_dofs(const size_t iT) const
{
    return m_high_cell_basis.at(iT)->dimension();
}

size_t XVEMCore::local_low_cell_dofs(const size_t iT) const
{
    return m_low_cell_basis.at(iT)->dimension();
}

size_t XVEMCore::local_high_edge_dofs(const size_t iE) const
{
    return m_high_edge_basis.at(iE)->dimension();
}

size_t XVEMCore::local_low_edge_dofs(const size_t iE) const
{
    return m_low_edge_basis.at(iE)->dimension();
    // return m_high_edge_basis.at(iE)->dimension() - 2;
}

size_t XVEMCore::local_dofs(const size_t iT) const
{
    size_t dofs = m_mesh->cell(iT)->n_vertices() + local_low_cell_dofs(iT);
    for (size_t iTF = 0; iTF < m_mesh->cell(iT)->n_edges(); ++iTF)
    {
        dofs += this->local_low_edge_dofs(m_mesh->cell(iT)->edge(iTF)->global_index());
    }
    return dofs;
}

size_t XVEMCore::global_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    size_t offset = 0;
    for (size_t i = 0; i < iT; ++i)
    {
        offset += this->local_low_cell_dofs(i);
    }
    return offset;
}

size_t XVEMCore::global_offset_E(const size_t iE) const
{
    assert(iE < m_mesh->n_edges());
    size_t offset = total_cell_dofs();
    for (size_t i = 0; i < iE; ++i)
    {
        offset += this->local_low_edge_dofs(i);
    }
    // offset += this->total_cell_dofs();
    return offset;
}

size_t XVEMCore::local_offset_E(const size_t iT, const size_t iTE) const
{
    assert(iTE < m_mesh->cell(iT)->n_edges());
    size_t offset = local_low_cell_dofs(iT);
    for (size_t i = 0; i < iTE; ++i)
    {
        size_t iE = m_mesh->cell(iT)->edge(i)->global_index();
        offset += this->local_low_edge_dofs(iE);
    }
    return offset;
}

size_t XVEMCore::local_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    return 0;
}

size_t XVEMCore::local_offset_V(const size_t iT, const size_t iTV) const
{
    assert((iT < m_mesh->n_cells()) && (iTV < m_mesh->cell(iT)->n_vertices()));
    return local_dofs(iT) - m_mesh->cell(iT)->n_vertices() + iTV;
}

size_t XVEMCore::global_offset_V(const size_t iV) const
{
    assert(iV < m_mesh->n_vertices());
    return total_cell_dofs() + total_edge_dofs() + iV;
}

size_t XVEMCore::total_cell_dofs() const
{
    size_t dofs = 0;
    for (size_t iT = 0; iT < m_mesh->n_cells(); ++iT)
    {
        dofs += this->local_low_cell_dofs(iT);
    }
    return dofs;
}
size_t XVEMCore::total_edge_dofs() const
{
    size_t dofs = 0;
    for (size_t iE = 0; iE < m_mesh->n_edges(); ++iE)
    {
        dofs += this->local_low_edge_dofs(iE);
    }
    return dofs;
}
size_t XVEMCore::total_dofs() const
{
    return total_cell_dofs() + total_edge_dofs() + m_mesh->n_vertices();
}
