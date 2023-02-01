#include <HybridCore.hpp>
#include <parallel_for.hpp>

using Quadrature::QuadratureRule;
using Quadrature::GaussLegendre1D;

HybridCore::HybridCore(
    MeshType *mesh_ptr, ///< A pointer to the loaded mesh
    size_t cell_deg,              ///< The degree of the cell polynomials
    size_t edge_deg,              ///< The degree of the edge polynomials
    bool use_threads,             ///< Optional argument to indicate if threads should be used
    bool ortho)

    : m_mesh(mesh_ptr),
      m_cell_deg(cell_deg),
      m_edge_deg(edge_deg),
      m_use_threads(use_threads),
      m_ortho(ortho),
      m_highorder_basis(mesh_ptr->n_cells()),
      m_cell_basis(mesh_ptr->n_cells()),
      m_edge_basis(mesh_ptr->n_edges())
{

    GaussLegendre1D low_quad(2 * (m_edge_deg + 1));
    GaussLegendre1D high_quad(30);

    for (auto &e : m_mesh->get_edges())
    {
        if(e->is_straight())
        {
            edge_quads.push_back(generate_quadrature_rule(e, low_quad));
        }
        else
        {
            edge_quads.push_back(generate_quadrature_rule(e, high_quad));
        }
    }

    for (auto &c : m_mesh->get_cells())
    {
        // std::vector<QuadratureRule<double>> local_edge_quads(c->n_edges());
        std::vector<QuadratureRule<double>> local_edge_quads;
        for (size_t iTF = 0; iTF < c->n_edges(); ++iTF)
        {
            // local_edge_quads[iTF] = edge_quads[c->edge(iTF)->global_index()];
            local_edge_quads.push_back(edge_quads[c->edge(iTF)->global_index()]);
        }
        cell_quads.push_back(generate_quadrature_rule(c, local_edge_quads, low_quad));
    }

    // boost::timer::cpu_timer construction_timer;
    // construction_timer.start();

    // std::cout << "\n[HybridCore] Constructing bases\n";

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

    // construction_timer.stop();
    // std::cout << "     Construction time = " << construction_timer.elapsed().wall * 1E-9 << "s\n";
}

HybridCore::CellBasisType HybridCore::construct_cell_basis(const size_t iT, const size_t degree)
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

    PolyMesh2D::Functional::VectorFunction2D transform(transform_val, transform_deriv);

    // Basis of monomials
    PolyMesh2D::Functional::ScalarBasis2D mono_basis = PolyMesh2D::Functional::MonomialScalarBasis(transform, degree);

    Eigen::MatrixXd basis_integrals = Eigen::MatrixXd::Zero(mono_basis.dimension(), mono_basis.dimension());

    for (size_t i = 0; i < mono_basis.dimension(); ++i)
    {
        for (size_t j = i; j < mono_basis.dimension(); ++j)
        {
            for (size_t iqn = 0; iqn < cell_quads[iT].size(); ++iqn)
            {
                basis_integrals(i, j) += cell_quads[iT][iqn].w * mono_basis.value(i, cell_quads[iT][iqn].x) * mono_basis.value(j, cell_quads[iT][iqn].x);
            }
            basis_integrals(j, i) = basis_integrals(i, j);
        }
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension());

    if (m_ortho)
    {
        double norm = std::sqrt(basis_integrals(0, 0));
        B(0, 0) = (1.0 / norm);
        basis_integrals.row(0) /= norm;
        Eigen::VectorXd tmp = basis_integrals.row(0).transpose();
        basis_integrals.col(0) = tmp;

        for (size_t i = 1; i < mono_basis.dimension(); ++i)
        {
            Eigen::RowVectorXd coeffs = Eigen::RowVectorXd::Zero(i);
            for (size_t j = 0; j < i; j++)
            {
                coeffs(j) = -basis_integrals(i, j);
            }

            for (size_t j = 0; j < i; j++)
            {
                basis_integrals.row(i) += coeffs(j) * basis_integrals.row(j);
            }

            norm = std::sqrt(basis_integrals(i, i));
            coeffs /= norm;
            basis_integrals.row(i) /= norm;

            tmp = basis_integrals.row(i).transpose();
            basis_integrals.col(i) = tmp;

            B.row(i).head(i) = coeffs * B.topLeftCorner(i, i);
            B(i, i) = (1.0 / norm);
        }
    }

    return CellBasisType(mono_basis, B);
}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

HybridCore::EdgeBasisType HybridCore::construct_edge_basis(const size_t iE)
{
    PolyMesh2D::Functional::ScalarBasis1D mono_basis;
    if (m_mesh->edge(iE)->is_straight())
    {
        mono_basis = PolyMesh2D::Functional::MonomialScalarBasis(m_mesh->edge(iE)->parameterisation(), m_edge_deg);
    }
    else
    {
        mono_basis = PolyMesh2D::Functional::CurvedEdgeBasis(m_mesh->edge(iE)->parameterisation(), m_edge_deg);
    }

    Eigen::MatrixXd basis_integrals = Eigen::MatrixXd::Zero(mono_basis.dimension(), mono_basis.dimension());

    for (size_t i = 0; i < mono_basis.dimension(); ++i)
    {
        for (size_t j = i; j < mono_basis.dimension(); ++j)
        {
            for (size_t iqn = 0; iqn < edge_quads[iE].size(); ++iqn)
            {
                basis_integrals(i, j) += edge_quads[iE][iqn].w * mono_basis.value(i, edge_quads[iE][iqn].x) * mono_basis.value(j, edge_quads[iE][iqn].x);
            }
            basis_integrals(j, i) = basis_integrals(i, j);
        }
    }

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
            if(to_remove)
            {
                mono_basis.remove_basis_function(basis_col - shift);
                ++shift;
            }
        }
        // std::cout << mono_basis.dimension() << "\n";

        new_basis_integrals = Eigen::MatrixXd::Zero(mono_basis.dimension(), mono_basis.dimension());

        for (size_t i = 0; i < mono_basis.dimension(); ++i)
        {
            for (size_t j = i; j < mono_basis.dimension(); ++j)
            {
                for (size_t iqn = 0; iqn < edge_quads[iE].size(); ++iqn)
                {
                    new_basis_integrals(i, j) += edge_quads[iE][iqn].w * mono_basis.value(i, edge_quads[iE][iqn].x) * mono_basis.value(j, edge_quads[iE][iqn].x);
                }
                new_basis_integrals(j, i) = new_basis_integrals(i, j);
            }
        }
    }
    else
    {
        new_basis_integrals = basis_integrals;
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Identity(mono_basis.dimension(), mono_basis.dimension());

    if (m_ortho)
    {
        double norm = std::sqrt(new_basis_integrals(0, 0));
        B(0, 0) = (1.0 / norm);
        new_basis_integrals.row(0) /= norm;
        Eigen::VectorXd tmp = new_basis_integrals.row(0).transpose();
        new_basis_integrals.col(0) = tmp;

        for (size_t i = 1; i < mono_basis.dimension(); ++i)
        {
            Eigen::RowVectorXd coeffs = Eigen::RowVectorXd::Zero(i);
            for (size_t j = 0; j < i; j++)
            {
                coeffs(j) = -new_basis_integrals(i, j);
            }

            for (size_t j = 0; j < i; j++)
            {
                new_basis_integrals.row(i) += coeffs(j) * new_basis_integrals.row(j);
            }

            norm = std::sqrt(new_basis_integrals(i, i));
            coeffs /= norm;
            new_basis_integrals.row(i) /= norm;

            tmp = new_basis_integrals.row(i).transpose();
            new_basis_integrals.col(i) = tmp;

            B.row(i).head(i) = coeffs * B.topLeftCorner(i, i);
            B(i, i) = (1.0 / norm);
        }
    }

    return EdgeBasisType(mono_basis, B);
}

Eigen::VectorXd HybridCore::restr(const Eigen::VectorXd &UVec, const size_t iT) const
{
    assert((size_t)UVec.size() == total_cell_dofs() + total_edge_dofs());

    Eigen::VectorXd UVec_T = Eigen::VectorXd::Zero(local_cell_dofs(iT) + local_boundary_dofs(iT));

    UVec_T.head(local_cell_dofs(iT)) = UVec.segment(global_offset_T(iT), local_cell_dofs(iT));

    for (size_t iTE = 0; iTE < m_mesh->cell(iT)->n_edges(); iTE++)
    {
        const size_t iE = m_mesh->cell(iT)->edge(iTE)->global_index();
        UVec_T.segment(local_cell_dofs(iT) + local_offset_E(iT, iTE), local_edge_dofs(iE)) = UVec.segment(total_cell_dofs() + global_offset_E(iE), local_edge_dofs(iE));
    }

    return UVec_T;
}

Eigen::VectorXd HybridCore::interpolate(const std::function<double(Eigen::Vector2d)> &func) const
{
    Eigen::VectorXd Ihk(Eigen::VectorXd::Zero(total_cell_dofs() + total_edge_dofs()));

    // Create cell bases
    std::function<void(size_t, size_t)> interpolate_cell_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            const size_t n_local_cell_dofs = local_cell_dofs(iT);

            Eigen::MatrixXd MTT = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_cell_dofs);
            Eigen::VectorXd u_integrated_against_basis_T = Eigen::VectorXd::Zero(n_local_cell_dofs);

            QuadratureRule<Eigen::Vector2d> quadT = cell_quads[iT];

            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                for (size_t i = 0; i < n_local_cell_dofs; ++i)
                {
                    // double basis_i_weight = quadT[iqn].w * m_cell_basis[iT]->value(i, quadT[iqn].x);
                    double basis_i_weight = quadT[iqn].w * m_highorder_basis[iT]->value(i, quadT[iqn].x);
                    u_integrated_against_basis_T(i) += basis_i_weight * func(quadT[iqn].x);
                    for (size_t j = 0; j < n_local_cell_dofs; ++j)
                    {
                        // MTT(i, j) += basis_i_weight * m_cell_basis[iT]->value(j, quadT[iqn].x);
                        MTT(i, j) += basis_i_weight * m_highorder_basis[iT]->value(j, quadT[iqn].x);
                    }
                }
            }

            Ihk.segment(global_offset_T(iT), n_local_cell_dofs) = MTT.inverse() * u_integrated_against_basis_T;
            // Ihk.segment(global_offset_T(iT), n_local_cell_dofs) = u_integrated_against_basis_T;
        }
    };
    parallel_for(m_mesh->n_cells(), interpolate_cell_terms, m_use_threads);

    // Create cell bases
    std::function<void(size_t, size_t)> interpolate_edge_terms = [&](size_t start, size_t end) -> void
    {
        for (size_t iF = start; iF < end; ++iF)
        {
            const size_t n_local_edge_dofs = local_edge_dofs(iF);

            Eigen::MatrixXd MFF = Eigen::MatrixXd::Zero(n_local_edge_dofs, n_local_edge_dofs);
            Eigen::VectorXd u_integrated_against_basis_F = Eigen::VectorXd::Zero(n_local_edge_dofs);

            QuadratureRule<double> quadF = edge_quads[iF];

            for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
            {
                double func_weight = quadF[iqn].w * func(m_mesh->edge(iF)->parameterisation().value(quadF[iqn].x));
                for (size_t i = 0; i < n_local_edge_dofs; ++i)
                {
                    u_integrated_against_basis_F(i) += func_weight * m_edge_basis[iF]->value(i, quadF[iqn].x);
                    for (size_t j = 0; j < n_local_edge_dofs; ++j)
                    {
                        MFF(i, j) += quadF[iqn].w * m_edge_basis[iF]->value(i, quadF[iqn].x) * m_edge_basis[iF]->value(j, quadF[iqn].x);
                    }
                }
            }

            Ihk.segment(total_cell_dofs() + global_offset_E(iF), n_local_edge_dofs) = MFF.inverse() * u_integrated_against_basis_F;
            // Ihk.segment(total_cell_dofs() + global_offset_E(iF), n_local_edge_dofs) = u_integrated_against_basis_F;
        }
    };
    parallel_for(m_mesh->n_edges(), interpolate_edge_terms, m_use_threads);

    return Ihk;
}

QuadratureRule<Eigen::Vector2d> HybridCore::quadT(const size_t iT) const
{
    assert(iT < cell_quads.size());
    return cell_quads[iT];
}
QuadratureRule<double> HybridCore::quadE(const size_t iE) const
{
    assert(iE < edge_quads.size());
    return edge_quads[iE];
}

HybridCore::MeshType *HybridCore::get_mesh() const
{
    return m_mesh;
}
size_t HybridCore::cell_degree() const
{
    return m_cell_deg;
}
size_t HybridCore::edge_degree() const
{
    return m_edge_deg;
}

HybridCore::CellBasisType *HybridCore::highorder_basis(const size_t iT) const
{
    assert(iT < m_highorder_basis.size());
    return m_highorder_basis[iT].get();
}
// HybridCore::CellBasisType *HybridCore::cell_basis(const size_t iT) const
// {
//     assert(iT < m_cell_basis.size());
//     return m_cell_basis[iT].get();
// }
HybridCore::EdgeBasisType *HybridCore::edge_basis(const size_t iE) const
{
    assert(iE < m_edge_basis.size());
    return m_edge_basis[iE].get();
}

size_t HybridCore::local_cell_dofs(const size_t iT) const
{
    // return m_cell_basis[iT]->dimension();
    // return m_cell_basis[iT]->dimension();
    return (m_cell_deg + 1) * (m_cell_deg + 2) / 2;
}
size_t HybridCore::local_highorder_dofs(const size_t iT) const
{
    return m_highorder_basis[iT]->dimension();
}
size_t HybridCore::local_edge_dofs(const size_t iE) const
{
    return m_edge_basis[iE]->dimension();
}

size_t HybridCore::local_boundary_dofs(const size_t iT) const
{
    size_t dofs = 0;
    for (size_t iTF = 0; iTF < m_mesh->cell(iT)->n_edges(); ++iTF)
    {
        dofs += this->local_edge_dofs(m_mesh->cell(iT)->edge(iTF)->global_index());
    }
    return dofs;
}

size_t HybridCore::global_offset_T(const size_t iT) const
{
    assert(iT < m_mesh->n_cells());
    size_t offset = 0;
    for (size_t i = 0; i < iT; ++i)
    {
        offset += this->local_cell_dofs(i);
    }
    return offset;
}
size_t HybridCore::global_offset_E(const size_t iE) const
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
size_t HybridCore::local_offset_E(const size_t iT, const size_t iTE) const
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
size_t HybridCore::local_offset_T(const size_t iT) const
{
    return 0;
}

size_t HybridCore::total_cell_dofs() const
{
    size_t dofs = 0;
    for (size_t iT = 0; iT < m_mesh->n_cells(); ++iT)
    {
        dofs += this->local_cell_dofs(iT);
    }
    return dofs;
}
size_t HybridCore::total_edge_dofs() const
{
    size_t dofs = 0;
    for (size_t iE = 0; iE < m_mesh->n_edges(); ++iE)
    {
        dofs += this->local_edge_dofs(iE);
    }
    return dofs;
}
