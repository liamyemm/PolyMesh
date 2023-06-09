#include "function.hpp"
#include "basis.hpp"
#include "common.hpp"
#include "QuadratureRule.hpp"
#include "generate_quadrature_rule.hpp"
#include "Mesh.hpp"

#include <iostream>

#include <Eigen/Dense>

#ifndef _QUADHANDLER_HPP
#define _QUADHANDLER_HPP

namespace PolyMesh2D
{
    /**
     * @brief Namespace for quadrature-related functionality.
     */
    namespace Quadrature
    {
        /**
         * @brief Class template for handling quadrature on 2D polygonal meshes.
         *
         * This class template provides functionality for integrating functions over
         * cells and edges of a polygonal mesh using Gauss-Legendre quadrature rules.
         * It works with both straight and curved polygonal meshes.
         *
         * @tparam MeshType The type of polygonal mesh to integrate over. Either StraightMesh::Mesh or CurvedMesh::Mesh.
         */
        template <typename MeshType>
        class QuadHandler
        {
        public:
            /**
             * @brief Constructs a new QuadHandler object.
             *
             * @param mesh      Pointer to the mesh object to integrate over.
             * @param cell_doe  Degree of exactness to use for cell integration.
             * @param edge_doe  Degree of exactness to use for edge integration.
             */
            QuadHandler(MeshType *mesh, size_t edge_doe, size_t cell_doe);

            /**
             * @brief Integrates a scalar-valued function over a cell.
             *
             * @param func  The function to integrate.
             * @param cell  Pointer to the cell to integrate over.
             * @return      The result of the integration.
             */
            double integrate(const Functional::Function<2, 1> &func, typename MeshType::CellType *cell) const;

            /**
             * @brief Integrates a scalar-valued function over an edge.
             *
             * @param func  The function to integrate.
             * @param edge  Pointer to the edge to integrate over.
             * @return      The result of the integration.
             */
            double integrate(const Functional::Function<1, 1> &func, typename MeshType::EdgeType *edge) const;

            /**
             * @brief Integrates a vector-valued function over a cell.
             *
             * @param func  The function to integrate.
             * @param cell  Pointer to the cell to integrate over.
             * @return      The result of the integration.
             */
            Eigen::Vector2d integrate(const Functional::Function<2, 2> &func, typename MeshType::CellType *cell) const;

            /**
             * @brief Calculates the L2 gram matrix of two bases over a given cell and returns the result as an Eigen matrix.
             *
             * This function computes the integral of the scalar product between of all pairs of basis functions in the bases `basis_1` and `basis_2`,
             * over a specified cell of the mesh. The output dimension of the basis functions should be specified as a template parameter.
             *
             * @tparam output_dim The output dimension of the basis functions.
             * @param basis_1 The first basis of type Functional::Basis<2, output_dim>.
             * @param basis_2 The second basis of type Functional::Basis<2, output_dim>.
             * @param cell Pointer to the cell of type MeshType::CellType over which the integration is performed.
             * @param sym Optional parameter to indicate whether or not the matrix is symmetric. Set to false by default.
             * @return The integrated scalar product matrix as an Eigen::MatrixXd.
             *
             * The function constructs an Eigen matrix, `mat`, with dimensions corresponding to the dimensions of
             * `basis_1` and `basis_2`. If the matrix is symmetric (i.e., `basis_1` is equal to `basis_2`),
             * the lower triangular values of the matrix are set to the corresponding upper triangular values
             * to avoid recomputing integrals. The integral of the scalar product between each pair of basis functions
             * is computed and assigned to the corresponding matrix elements.
             */
            template <unsigned output_dim>
            Eigen::MatrixXd l2_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::CellType *cell, bool sym = false) const;

            template <unsigned output_dim>
            Eigen::MatrixXd l2_product(const Functional::Basis<1, output_dim> &basis_1, const Functional::Basis<1, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym = false) const;

            template <unsigned output_dim>
            Eigen::MatrixXd l2_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym = false) const;

            template <unsigned output_dim>
            Eigen::MatrixXd l2_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<1, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym = false) const;

            template <unsigned output_dim>
            Eigen::MatrixXd l2_product(const Functional::Basis<1, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym = false) const;

            template <typename BasisType1, typename BasisType2, typename DomainType>
            Eigen::MatrixXd l2_product(const Functional::Family<BasisType1> &family_1, const Functional::Family<BasisType2> &family_2, DomainType *cell_or_edge, bool sym = false) const;

            template <unsigned input_dim, unsigned output_dim, typename DomainType>
            Eigen::VectorXd l2_product(const Functional::Function<input_dim, output_dim> &func, const Functional::Basis<input_dim, output_dim> &basis, DomainType *cell_or_edge) const;

            template <unsigned input_dim, unsigned output_dim, typename DomainType>
            Eigen::VectorXd l2_product(const Functional::Function<input_dim, output_dim> &func, const Functional::Family<Functional::Basis<input_dim, output_dim>> &family, DomainType *cell_or_edge) const;

            template <unsigned output_dim>
            Eigen::MatrixXd h1_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::CellType *cell, bool sym = false) const;

            template <typename BasisType1, typename BasisType2>
            Eigen::MatrixXd h1_product(const Functional::Family<BasisType1> &family_1, const Functional::Family<BasisType2> &family_2, typename MeshType::CellType *cell, bool sym = false) const;

            const QuadratureRule<Eigen::Vector2d> &get_cell_quad(size_t index) const
            {
                return cell_quads.at(index);
            }

        private:
            std::vector<QuadratureRule<double>> edge_quads;          ///< Vector of Gauss-Legendre quadrature rules for edge integration.
            std::vector<QuadratureRule<Eigen::Vector2d>> cell_quads; ///< Vector of Gauss-Legendre quadrature rules for cell integration.

            const size_t _edge_doe; ///< The degree of exactness to use for edge integration.
            const size_t _cell_doe; ///< The degree of exactness to use for cell integration.
        };

        template <typename MeshType>
        template <unsigned output_dim>
        Eigen::MatrixXd QuadHandler<MeshType>::l2_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::CellType *cell, bool sym) const
        {
            Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(basis_1.dimension(), basis_2.dimension());

            if (sym) // if symmetric, set lower triangular values to upper triangular values to avoid recomputing integrals
            {
                // std::cout << "\n\nhell1\n\n";
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    for (size_t j = i; j < basis_2.dimension(); ++j)
                    {
                        mat(i, j) = this->integrate(Functional::scalar_product(basis_1.function(i), basis_2.function(j)), cell);
                        mat(j, i) = mat(i, j);
                    }
                }
            }
            else
            {
                // std::cout << "\n\nhellO\n\n";
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    for (size_t j = 0; j < basis_2.dimension(); ++j)
                    {
                        mat(i, j) = this->integrate(Functional::scalar_product(basis_1.function(i), basis_2.function(j)), cell);
                    }
                }
            }

            return mat;
        }

        template <typename MeshType>
        template <unsigned input_dim, unsigned output_dim, typename DomainType>
        Eigen::VectorXd QuadHandler<MeshType>::l2_product(const Functional::Function<input_dim, output_dim> &func, const Functional::Basis<input_dim, output_dim> &basis, DomainType *cell_or_edge) const
        {
            Eigen::VectorXd vec(basis.dimension());

            for (size_t i = 0; i < basis.dimension(); ++i)
            {
                vec(i) = this->integrate(Functional::scalar_product(func, basis.function(i)), cell_or_edge);
            }
            return vec;
        }

        template <typename MeshType>
        template <unsigned input_dim, unsigned output_dim, typename DomainType>
        Eigen::VectorXd QuadHandler<MeshType>::l2_product(const Functional::Function<input_dim, output_dim> &func, const Functional::Family<Functional::Basis<input_dim, output_dim>> &family, DomainType *cell_or_edge) const
        {
            return family.matrix() * this->l2_product(func, family.ancestor(), cell_or_edge);
        }

        template <typename MeshType>
        template <typename BasisType1, typename BasisType2, typename DomainType>
        Eigen::MatrixXd QuadHandler<MeshType>::l2_product(const Functional::Family<BasisType1> &family_1, const Functional::Family<BasisType2> &family_2, DomainType *cell_or_edge, bool sym) const
        {
            Eigen::MatrixXd basis_gram_mat = this->l2_product(family_1.ancestor(), family_2.ancestor(), cell_or_edge, sym);

            return family_1.matrix() * basis_gram_mat * family_2.matrix().transpose();
        }

        template <typename MeshType>
        template <typename BasisType1, typename BasisType2>
        Eigen::MatrixXd QuadHandler<MeshType>::h1_product(const Functional::Family<BasisType1> &family_1, const Functional::Family<BasisType2> &family_2, typename MeshType::CellType *cell, bool sym) const
        {
            Eigen::MatrixXd basis_gram_mat = this->h1_product(family_1.ancestor(), family_2.ancestor(), cell, sym);

            return family_1.matrix() * basis_gram_mat * family_2.matrix().transpose();
        }

        template <typename MeshType>
        template <unsigned output_dim>
        Eigen::MatrixXd QuadHandler<MeshType>::h1_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::CellType *cell, bool sym) const
        {
            Eigen::MatrixXd mat(basis_1.dimension(), basis_2.dimension());

            if (sym) // if symmetric, set lower triangular values to upper triangular values to avoid recomputing integrals
            {
                // std::cout << "\n\nIamsim\n\n";
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    for (size_t j = i; j < basis_2.dimension(); ++j)
                    {
                        mat(j, i) = mat(i, j) = this->integrate(Functional::gradient_product(basis_1.function(i), basis_2.function(j)), cell);
                    }
                }
            }
            else
            {
                // std::cout << "\n\nIamnotsim\n\n";
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    for (size_t j = 0; j < basis_2.dimension(); ++j)
                    {
                        mat(i, j) = this->integrate(Functional::gradient_product(basis_1.function(i), basis_2.function(j)), cell);
                    }
                }
            }

            return mat;
        }

        template <typename MeshType>
        template <unsigned output_dim>
        Eigen::MatrixXd QuadHandler<MeshType>::l2_product(const Functional::Basis<1, output_dim> &basis_1, const Functional::Basis<1, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym) const
        {
            Eigen::MatrixXd mat(basis_1.dimension(), basis_2.dimension());

            if (sym) // if symmetric, set lower triangular values to upper triangular values to avoid recomputing integrals
            {
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    for (size_t j = i; j < basis_2.dimension(); ++j)
                    {
                        mat(j, i) = mat(i, j) = this->integrate(Functional::scalar_product(basis_1.function(i), basis_2.function(j)), edge);
                    }
                }
            }
            else
            {
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    for (size_t j = 0; j < basis_2.dimension(); ++j)
                    {
                        mat(i, j) = this->integrate(Functional::scalar_product(basis_1.function(i), basis_2.function(j)), edge);
                    }
                }
            }

            return mat;
        }

        template <typename MeshType>
        template <unsigned output_dim>
        Eigen::MatrixXd QuadHandler<MeshType>::l2_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym) const
        {
            Eigen::MatrixXd mat(basis_1.dimension(), basis_2.dimension());

            if (sym) // if symmetric, set lower triangular values to upper triangular values to avoid recomputing integrals
            {
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    auto func_i_trace = Functional::trace(basis_1.function(i), edge);
                    for (size_t j = i; j < basis_2.dimension(); ++j)
                    {
                        auto func_j_trace = Functional::trace(basis_2.function(j), edge);
                        mat(j, i) = mat(i, j) = this->integrate(Functional::scalar_product(func_i_trace, func_j_trace), edge);
                    }
                }
            }
            else
            {
                for (size_t i = 0; i < basis_1.dimension(); ++i)
                {
                    auto func_i_trace = Functional::trace(basis_1.function(i), edge);
                    for (size_t j = 0; j < basis_2.dimension(); ++j)
                    {
                        auto func_j_trace = Functional::trace(basis_2.function(j), edge);
                        mat(i, j) = this->integrate(Functional::scalar_product(func_i_trace, func_j_trace), edge);
                    }
                }
            }

            return mat;
        }

        template <typename MeshType>
        template <unsigned output_dim>
        Eigen::MatrixXd QuadHandler<MeshType>::l2_product(const Functional::Basis<2, output_dim> &basis_1, const Functional::Basis<1, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym) const
        {
            if (sym == true)
            {
                throw "Error - this gram matrix cannot be symmetric\n";
            }
            Eigen::MatrixXd mat(basis_1.dimension(), basis_2.dimension());

            for (size_t i = 0; i < basis_1.dimension(); ++i)
            {
                auto func_i_trace = Functional::trace(basis_1.function(i), edge->parameterisation());
                for (size_t j = 0; j < basis_2.dimension(); ++j)
                {
                    mat(i, j) = this->integrate(Functional::scalar_product(func_i_trace, basis_2.function(j)), edge);
                }
            }

            return mat;
        }

        template <typename MeshType>
        template <unsigned output_dim>
        Eigen::MatrixXd QuadHandler<MeshType>::l2_product(const Functional::Basis<1, output_dim> &basis_1, const Functional::Basis<2, output_dim> &basis_2, typename MeshType::EdgeType *edge, bool sym) const
        {
            return this->l2_product(basis_2, basis_1, edge, sym).transpose();
        }
    }
}

#endif