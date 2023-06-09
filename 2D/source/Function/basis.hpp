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

        /**
         * @brief Templated class for defining a basis of functions of arbitrary input and output dimensions
         *
         * @tparam input_dim Number of input dimensions
         * @tparam output_dim Number of output dimensions
         */
        template <unsigned input_dim, unsigned output_dim>
        class Basis
        {
        public:
            /** Alias for the input type of basis functions **/
            typedef typename Function<input_dim, output_dim>::InputType InputType;

            /** Alias for the output type of basis functions **/
            typedef typename Function<input_dim, output_dim>::OutputType OutputType;

            /** Alias for the derivative type of basis functions **/
            typedef typename Function<input_dim, output_dim>::DerivativeType DerivativeType;

            /**
             * Number of input dimensions
             **/
            static const unsigned input_dimension = input_dim;

            /**
             * Number of output dimensions
             **/
            static const unsigned output_dimension = output_dim;

            /**
             * Constructor
             **/
            Basis() {}

            /**
             * Destructor
             **/
            ~Basis() {}

            /**
             * Constructor with a set of basis functions
             * @param bases A vector of basis functions
             **/
            Basis(const std::vector<Function<input_dim, output_dim>> &bases) : __basis(bases) {}

            /**
             * Get the number of basis functions
             * @return size_t Number of basis functions
             **/
            size_t dimension() const
            {
                return __basis.size();
            }

            /**
             * Add a basis function
             * @param func The basis function to add
             **/
            void add_basis_function(const Function<input_dim, output_dim> &func)
            {
                __basis.push_back(func);
            }

            /**
             * Remove a basis function at a given index
             * @param i The index of the basis function to remove
             **/
            void remove_basis_function(const size_t i) // remove basis function at ith index
            {
                assert(i < dimension());
                __basis.erase(std::begin(__basis) + i);
            }

            /**
             * Evaluate the i-th basis function at point x
             * @param i The index of the basis function to evaluate
             * @param x The input point at which to evaluate the basis function
             * @return OutputType The output value of the basis function
             **/
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

            /**
             * Evaluate the derivative of the i-th basis function at point x
             * @param i The index of the basis function to evaluate
             * @param x The input point at which to evaluate the basis function
             * @return DerivativeType The derivative value of the basis function
             **/
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

            /**
             * Returns a const reference to the i-th basis Function.
             * @param i The index of the basis Function to retrieve.
             * @return A const reference to the i-th basis Function.
             * @throws An exception if the index is out of bounds.
             * This function returns a const reference to the i-th basis Function in this basis set. If the index is out of bounds (i.e. greater than or equal to the dimension of this basis set), an exception is thrown with an error message indicating the index that was out of bounds.
             **/
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

        /**
         * @brief A class representing a family of functions expressed in a given basis with coefficients stored in a matrix. This class allows the evaluation of functions in the family and their derivatives at any given point in the domain of the functions. The family is represented by a BasisType object which defines the basis in which the family is expressed, and a coefficient matrix whose i-th line contains the coefficient of the expansion of the i-th function of the family in the basis.
         * @tparam BasisType Type of the basis object used to represent the family
         **/
        template <typename BasisType>
        class Family
        {
        public:
            /** An alias for the input type of the basis **/
            typedef typename BasisType::InputType InputType;

            /** An alias for the output type of the basis **/
            typedef typename BasisType::OutputType OutputType;

            /** An alias for the derivative type of the basis **/
            typedef typename BasisType::DerivativeType DerivativeType; 

            Family( ) {}

            /**
             * Constructor for the Family class
             * @param basis The basis in which the family is expressed
             * @param matrix The coefficient matrix whose i-th line contains the coefficient of the expansion of the i-th Function of the family in the basis
             **/
            Family(
                const BasisType &basis,       
                const Eigen::MatrixXd &matrix 
                )
                : m_basis(basis),
                  m_matrix(matrix)
            {
                assert((size_t)m_matrix.cols() == m_basis.dimension());
                assert((size_t)m_matrix.rows() == m_basis.dimension());
            }

            /**
             * Returns the dimension of the family
             * This is actually the number of Functions in the family, not necessarily linearly independent
             * @return size_t The number of functions in the family
             **/
            inline size_t dimension() const
            {
                return m_basis.dimension();
            }

            /**
             * Adds a basis function to the family
             * @param func The basis function to add
             **/
            void add_basis_function(const Function<BasisType::input_dimension, BasisType::output_dimension> &func)
            {
                m_basis.add_basis_function(func);
            }

            /**
             * Removes a basis function from the family
             * @param i The index of the basis function to remove
             **/
            void remove_basis_function(const size_t i) // remove basis function at ith index
            {
                assert(i < dimension());
                m_basis.remove_basis_function(i);
            }

            /**
             * Resets the coefficient matrix of the family
             * @param matrix The new coefficient matrix
             **/
            void reset_matrix(const Eigen::MatrixXd &matrix)
            {
                assert((size_t)matrix.cols() == m_basis.dimension());
                m_matrix = matrix;
            }

            /**
             * Evaluate the i-th function in the family at a given point x
             * @param i The index of the function to evaluate
             * @param x The point at which to evaluate the function
             * @return OutputType The value of the function at the given point
             **/
            OutputType value(const size_t i, const InputType &x) const
            {
                OutputType f = m_matrix(i, 0) * m_basis.value(0, x);
                for (auto j = 1; j < m_matrix.cols(); j++)
                {
                    f += m_matrix(i, j) * m_basis.value(j, x);
                } // for j
                return f;
            }

            /**
             * Evaluate the gradient of the i-th Function at point x.
             * @param i The index of the Function to evaluate.
             * @param x The point at which to evaluate the Function.
             * @return The gradient of the i-th Function at point x.
             **/
            DerivativeType derivative(const size_t i, const InputType &x) const
            {
                DerivativeType G = m_matrix(i, 0) * m_basis.derivative(0, x);
                for (auto j = 1; j < m_matrix.cols(); j++)
                {
                    G += m_matrix(i, j) * m_basis.derivative(j, x);
                } // for j
                return G;
            }

            /**
             * Return the coefficient matrix.
             * @return The coefficient matrix.
             **/
            const Eigen::MatrixXd &matrix() const
            {
                return m_matrix;
            }

            /**
             * Return the ancestor basis functions.
             * @return The ancestor basis functions.
             **/
            const BasisType &ancestor() const
            {
                return m_basis;
            }

        private:
            BasisType m_basis;
            Eigen::MatrixXd m_matrix;
        };

        /**
         * Alias for a 1D scalar Basis.
         **/
        using ScalarBasis1D = Basis<1, 1>;

        /**
         * Alias for a 2D scalar Basis.
         **/
        using ScalarBasis2D = Basis<2, 1>;

        /**
         * Alias for a 1D vector Basis.
         **/
        using VectorBasis1D = Basis<1, 2>;

        /**
         * Alias for a 2D vector Basis.
         **/
        using VectorBasis2D = Basis<2, 2>;

        /**
         * Alias for a 1D scalar Family.
         **/
        using ScalarFamily1D = Family<ScalarBasis1D>;

        /**
         * Alias for a 2D scalar Family.
         **/
        using ScalarFamily2D = Family<ScalarBasis2D>;

        /**
         * Alias for a 1D vector Family.
         */
        using VectorFamily1D = Family<VectorBasis1D>;

        /**
         * Alias for a 2D vector Family.
         **/
        using VectorFamily2D = Family<VectorBasis2D>;

        //@}
    }
}
#endif
