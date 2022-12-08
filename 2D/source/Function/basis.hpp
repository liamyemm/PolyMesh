#include <function.hpp>
#include <Eigen/Dense>

#ifndef _BASIS_HPP
#define _BASIS_HPP

namespace PolyMesh2D
{
    namespace Functional
    {
        /*!
         * \addtogroup Functional
         * @{
         */

        template <unsigned input_dim, unsigned output_dim>
        class Basis
        {
        public:
            typedef typename Function<input_dim, output_dim>::InputType InputType;
            typedef typename Function<input_dim, output_dim>::OutputType OutputType;
            typedef typename Function<input_dim, output_dim>::DerivativeType DerivativeType;

            Basis() {}
            ~Basis() {}
            Basis(const std::vector<Function<input_dim, output_dim>> &bases) : __basis(bases) {}

            size_t dimension() const
            {
                return __basis.size();
            }

            void add_basis_function(const Function<input_dim, output_dim> &func)
            {
                return __basis.push_back(func);
            }

            void remove_basis_function(const size_t i) // remove basis function at ith index
            {
                assert(i < dimension());
                __basis.erase(std::begin(__basis) + i);
            }

            OutputType value(const size_t i, const InputType &x) const
            {
                if (i < dimension())
                {
                    return __basis[i].value(x);
                }
                else
                {
                    throw "ERROR - index " + std::to_string(i) + " out of bounds\n";
                }
            }

            DerivativeType derivative(const size_t i, const InputType &x) const
            {
                if (i < dimension())
                {
                    return __basis[i].derivative(x);
                }
                else
                {
                    throw "ERROR - index " + std::to_string(i) + " out of bounds\n";
                }
            }

            const Function<input_dim, output_dim> &function(const size_t i) const
            {
                if (i < dimension())
                {
                    return __basis[i];
                }
                else
                {
                    throw "ERROR - index " + std::to_string(i) + " out of bounds\n";
                }
            }

        private:
            std::vector<Function<input_dim, output_dim>> __basis;
        };

        template <typename BasisType>
        class Family
        {
        public:
            typedef typename BasisType::InputType InputType;
            typedef typename BasisType::OutputType OutputType;
            typedef typename BasisType::DerivativeType DerivativeType;
            /// Constructor
            Family(
                const BasisType &basis,       ///< The basis in which the family is expressed
                const Eigen::MatrixXd &matrix ///< The coefficient matrix whose i-th line contains the coefficient of the expansion of the i-th Function of the family in the basis
                )
                : m_basis(basis),
                  m_matrix(matrix)
            {
                assert((size_t)m_matrix.cols() == m_basis.dimension());
                assert((size_t)m_matrix.rows() == m_basis.dimension());
            }

            /// Dimension of the family. This is actually the number of Functions in the family, not necessarily linearly independent
            inline size_t dimension() const
            {
                return m_basis.dimension();
            }

            /// Evaluate the i-th Function at point x
            OutputType value(const size_t i, const InputType &x) const
            {
                OutputType f = m_matrix(i, 0) * m_basis.value(0, x);
                for (auto j = 1; j < m_matrix.cols(); j++)
                {
                    f += m_matrix(i, j) * m_basis.value(j, x);
                } // for j
                return f;
            }

            /// Evaluate the gradient of the i-th Function at point x
            DerivativeType derivative(const size_t i, const InputType &x) const
            {
                DerivativeType G = m_matrix(i, 0) * m_basis.derivative(0, x);
                for (auto j = 1; j < m_matrix.cols(); j++)
                {
                    G += m_matrix(i, j) * m_basis.derivative(j, x);
                } // for j
                return G;
            }

            /// Return the coefficient matrix
            inline const Eigen::MatrixXd &matrix() const
            {
                return m_matrix;
            }

            /// Return the ancestor
            inline const BasisType &ancestor() const
            {
                return m_basis;
            }

        private:
            BasisType m_basis;
            Eigen::MatrixXd m_matrix;
        };

        using ScalarBasis1D = Basis<1, 1>;
        using ScalarBasis2D = Basis<2, 1>;
        using VectorBasis1D = Basis<1, 2>;
        using VectorBasis2D = Basis<2, 2>;

        using ScalarFamily1D = Family<ScalarBasis1D>;
        using ScalarFamily2D = Family<ScalarBasis2D>;
        using VectorFamily1D = Family<VectorBasis1D>;
        using VectorFamily2D = Family<VectorBasis2D>;

        //@}
    }
}
#endif
