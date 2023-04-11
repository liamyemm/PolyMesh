#include <Mesh.hpp>
#include <cmath>
#include <function.hpp>
#include <common.hpp>
#include <basis.hpp>
#include <QuadratureRule.hpp>
#include <memory>
#include <iostream>

#ifndef _HYBRIDCORE_HPP
#define _HYBRIDCORE_HPP

//using Quadrature::QuadratureRule;
//using Quadrature::GaussLegendre1D;

class HybridCore
{
public:

using CellBasisType = PolyMesh2D::Functional::ScalarFamily2D;
using EdgeBasisType = PolyMesh2D::Functional::ScalarFamily1D;

using MeshType = PolyMesh2D::CurvedMesh::Mesh;

    HybridCore(
        MeshType *mesh_ptr,  ///< A pointer to the loaded mesh
        size_t cell_deg, ///< The degree of the cell polynomials
        size_t edge_deg, ///< The degree of the edge polynomials
        bool use_threads = true, ///< Optional argument to indicate if threads should be used
        bool ortho = true);

    Eigen::VectorXd restr(const Eigen::VectorXd &UVec, const size_t iT) const;

    Eigen::VectorXd interpolate(const std::function<double(Eigen::Vector2d)> &func) const;

    MeshType *get_mesh() const;
    size_t cell_degree() const;
    size_t edge_degree() const;

    CellBasisType *highorder_basis(const size_t iT) const;
    CellBasisType *cell_basis(const size_t iT) const;
    EdgeBasisType *edge_basis(const size_t iE) const;

    void enrich_highorder_basis(const size_t iT, const PolyMesh2D::Functional::ScalarFunction2D &func);
    void enrich_cell_basis(const size_t iT, const PolyMesh2D::Functional::ScalarFunction2D &func);
    void enrich_edge_basis(const size_t iE, const PolyMesh2D::Functional::ScalarFunction1D &func);

    Quadrature::QuadratureRule<Eigen::Vector2d> quadT(const size_t iT) const;
    Quadrature::QuadratureRule<double> quadE(const size_t iE) const;

    size_t local_cell_dofs(const size_t iT) const;
    size_t local_highorder_dofs(const size_t iT) const;
    size_t local_edge_dofs(const size_t iE) const;
    size_t local_boundary_dofs(const size_t iT) const;

    size_t global_offset_T(const size_t iT) const;
    size_t global_offset_E(const size_t iE) const;
    size_t local_offset_E(const size_t iT, const size_t iTE) const;
    size_t local_offset_T(const size_t iT) const;

    size_t total_cell_dofs() const;
    size_t total_edge_dofs() const;

private:
    MeshType *m_mesh; // Pointer to mesh data
    size_t m_cell_deg;
    size_t m_edge_deg;

    bool m_use_threads;
    bool m_ortho;

    std::vector<std::unique_ptr<CellBasisType>> m_highorder_basis;
    std::vector<std::unique_ptr<CellBasisType>> m_cell_basis;
    std::vector<std::unique_ptr<EdgeBasisType>> m_edge_basis;

    std::vector<Quadrature::QuadratureRule<double>> edge_quads;
    std::vector<Quadrature::QuadratureRule<Eigen::Vector2d>> cell_quads;

    // Function2D::ScalarFamily2D construct_highorder_basis(const size_t iT);
    CellBasisType construct_cell_basis(const size_t iT, const size_t degree);
    EdgeBasisType construct_edge_basis(const size_t iE);
};

#endif
