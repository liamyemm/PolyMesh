#include <Mesh.hpp>
#include <cmath>
#include <function.hpp>
#include <common.hpp>
#include <basis.hpp>
#include <QuadratureRule.hpp>
#include <memory>
#include <iostream>

#include "QuadHandler.hpp"

#ifndef _XVEMCORE_HPP
#define _XVEMCORE_HPP

namespace PolyMesh2D
{
    class XVEMCore
    {
    public:
        using CellBasisType = Functional::ScalarFamily2D;
        using EdgeBasisType = Functional::ScalarFamily1D;

        using MeshType = CurvedMesh::Mesh;

        XVEMCore(
            MeshType *mesh_ptr,      ///< A pointer to the loaded mesh
            size_t l,         ///< The degree of the cell polynomials
            size_t k,         ///< The degree of the edge polynomials
            bool use_threads = true, ///< Optional argument to indicate if threads should be used
            bool ortho = true);

        Eigen::VectorXd restr(const Eigen::VectorXd &UVec, const size_t iT) const;

        Eigen::VectorXd interpolate(const Functional::Function<2, 1> &func) const;

        MeshType *get_mesh() const;

        size_t get_k() const;
        size_t get_l() const;

        CellBasisType *high_cell_basis(const size_t iT) const;
        CellBasisType *low_cell_basis(const size_t iT) const;

        EdgeBasisType *high_edge_basis(const size_t iE) const;
        EdgeBasisType *low_edge_basis(const size_t iE) const;

        void enrich_cell_basis(const size_t iT, const Functional::ScalarFunction2D &func);
        void enrich_edge_basis(const size_t iE, const Functional::ScalarFunction1D &func);

        Eigen::VectorXd l2_low_cell_projection(const size_t iT, const Functional::ScalarFunction2D &func) const;
        Eigen::VectorXd l2_low_edge_projection(const size_t iE, const Functional::ScalarFunction1D &func) const;

        Eigen::VectorXd elliptic_projection(const size_t iT, const Functional::ScalarFunction2D &func) const;

        size_t local_high_cell_dofs(const size_t iT) const;
        size_t local_low_cell_dofs(const size_t iT) const;

        size_t local_high_edge_dofs(const size_t iE) const;
        size_t local_low_edge_dofs(const size_t iE) const;

        size_t local_dofs(const size_t iT) const;

        size_t local_offset_E(const size_t iT, const size_t iTE) const;
        size_t local_offset_T(const size_t iT) const;
        size_t local_offset_V(const size_t iT, const size_t iTV) const;


        size_t global_offset_T(const size_t iT) const;
        size_t global_offset_E(const size_t iE) const;
        size_t global_offset_V(const size_t iV) const;

        size_t total_cell_dofs() const;
        size_t total_edge_dofs() const;

        size_t total_dofs() const;

        const Quadrature::QuadHandler<CurvedMesh::Mesh> &get_quad_handle() const;

        Eigen::MatrixXd continuous_reconstruction(const size_t iT, const size_t iTE) const;
        Eigen::MatrixXd potential_reconstruction(const size_t iT) const;

        Eigen::MatrixXd edge_restrict_mat(const size_t iT, const size_t iTE) const; // matrix which maps the local dofs of an element and its boundary to the dofs of an edge and its vertices

    private:
        MeshType *m_mesh; // Pointer to mesh data

        size_t m_l;
        size_t m_k;

        Quadrature::QuadHandler<CurvedMesh::Mesh> m_quad_handle;

        bool m_use_threads;
        bool m_ortho;

        std::vector<std::unique_ptr<CellBasisType>> m_high_cell_basis;
        std::vector<std::unique_ptr<CellBasisType>> m_low_cell_basis;

        std::vector<std::unique_ptr<EdgeBasisType>> m_high_edge_basis;
        std::vector<std::unique_ptr<EdgeBasisType>> m_low_edge_basis;

        // Function2D::ScalarFamily2D construct_highorder_basis(const size_t iT);
        CellBasisType construct_cell_basis(const size_t iT, const size_t degree) const;
        EdgeBasisType construct_edge_basis(const size_t iE, const int degree) const;
    };
}

#endif
