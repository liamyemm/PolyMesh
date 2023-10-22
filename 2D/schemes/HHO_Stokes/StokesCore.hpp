#include <Mesh.hpp>
#include <cmath>
#include <function.hpp>
#include <common.hpp>
#include <basis.hpp>
#include <QuadratureRule.hpp>
#include <memory>
#include <iostream>

#include "QuadHandler.hpp"

#ifndef _STOKESCORE_HPP
#define _STOKESCORE_HPP

namespace PolyMesh2D
{

    class StokesCore
    {
    public:
        using CellBasisType = Functional::VectorFamily2D;
        using EdgeBasisType = Functional::VectorFamily1D;

        using ScalarBasisType = Functional::ScalarFamily2D;

        using MeshType = CurvedMesh::Mesh;

        StokesCore(
            MeshType *mesh_ptr,      ///< A pointer to the loaded mesh
            size_t cell_deg,         ///< The degree of the cell polynomials
            size_t edge_deg,         ///< The degree of the edge polynomials
            bool use_threads = true, ///< Optional argument to indicate if threads should be used
            bool ortho = true);

        Eigen::VectorXd velocity_restr(const Eigen::VectorXd &UVec, const size_t iT) const;
        Eigen::VectorXd pressure_restr(const Eigen::VectorXd &UVec, const size_t iT) const;

        Eigen::VectorXd interpolate(const Functional::Function<2, 2> &u_func, const Functional::Function<2, 1> &p_func) const;

        MeshType *get_mesh() const;
        size_t cell_degree() const;
        size_t edge_degree() const;

        CellBasisType *highorder_basis(const size_t iT) const;
        CellBasisType *cell_basis(const size_t iT) const;
        EdgeBasisType *edge_basis(const size_t iE) const;

        Functional::Function<2, 1> bubble(const CurvedMesh::Cell *cell) const;

        Functional::Function<2, 2> vector_bubble(const CurvedMesh::Cell *cell) const;
        Functional::Function<2, 1> vector_bubble_divergence(const CurvedMesh::Cell *cell) const;

        ScalarBasisType *pressure_basis(const size_t iT) const;

        void enrich_highorder_basis(const size_t iT, const Functional::VectorFunction2D &func);
        void enrich_cell_basis(const size_t iT, const Functional::VectorFunction2D &func);
        void enrich_edge_basis(const size_t iE, const Functional::VectorFunction1D &func);

        void enrich_pressure_basis(const size_t iT, const Functional::ScalarFunction2D &func);

        Eigen::VectorXd pressure_robust_RHS(const size_t iT, const Functional::VectorFunction2D &source, const std::vector<Functional::Function<2, 1>> &RTN_enrichment = std::vector<Functional::Function<2, 1>>(), const std::vector<Functional::Function<2, 2>> &laplace_cell_enrichment = std::vector<Functional::Function<2, 2>>()) const;

        Eigen::VectorXd l2_cell_projection(const size_t iT, const Functional::VectorFunction2D &func) const;
        Eigen::VectorXd l2_edge_projection(const size_t iE, const Functional::VectorFunction1D &func) const;
        Eigen::VectorXd l2_highorder_projection(const size_t iT, const Functional::VectorFunction2D &func) const;
        Eigen::VectorXd l2_pressure_projection(const size_t iT, const Functional::ScalarFunction2D &func) const;

        Eigen::VectorXd elliptic_projection(const size_t iT, const Functional::VectorFunction2D &func) const;

        const Quadrature::QuadHandler<CurvedMesh::Mesh> &get_quad_handle() const;

        size_t local_cell_dofs(const size_t iT) const;
        size_t local_highorder_dofs(const size_t iT) const;
        size_t local_edge_dofs(const size_t iE) const;
        size_t local_boundary_dofs(const size_t iT) const;

        size_t local_pressure_dofs(const size_t iT) const;

        size_t global_offset_T(const size_t iT) const;
        size_t global_offset_E(const size_t iE) const;
        size_t local_offset_E(const size_t iT, const size_t iTE) const;
        size_t local_offset_T(const size_t iT) const;

        size_t global_pressure_offset_T(const size_t iT) const;
        size_t local_pressure_offset_T(const size_t iT) const;

        size_t total_cell_dofs() const;
        size_t total_edge_dofs() const;
        size_t total_pressure_dofs() const;

    private:
        MeshType *m_mesh; // Pointer to mesh data
        size_t m_cell_deg;
        size_t m_edge_deg;

        Quadrature::QuadHandler<CurvedMesh::Mesh> m_quad_handle;

        bool m_use_threads;
        bool m_ortho;

        std::vector<std::unique_ptr<CellBasisType>> m_highorder_basis;
        std::vector<std::unique_ptr<CellBasisType>> m_cell_basis;
        std::vector<std::unique_ptr<EdgeBasisType>> m_edge_basis;
        std::vector<std::unique_ptr<ScalarBasisType>> m_pressure_basis;

        // Function2D::ScalarFamily2D construct_highorder_basis(const size_t iT);
        CellBasisType construct_cell_basis(const size_t iT, const size_t degree) const;
        EdgeBasisType construct_edge_basis(const size_t iE) const;
        ScalarBasisType construct_pressure_basis(const size_t iT, const size_t degree) const;

        PolyMesh2D::Functional::ScalarFamily1D construct_scalar_edge_basis(const size_t iE, const size_t degree) const;
    };
}

#endif
