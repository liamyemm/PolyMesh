// Provides wrapper classes and methods for abstract functions that map R^n to R^m
//
// Author: Liam Yemm (liam.yemm@monash.edu)

#ifndef _FUNCTION_HPP
#define _FUNCTION_HPP

#include <functional> // std::function
#include <Eigen/Dense> // Eigen::Matrix

/*!
* @defgroup Functional
* @brief Classes and methods for generic functions and bases of function spaces.
*/

/** The global namespace for all classes and methods in the Polymesh2D project. **/
namespace PolyMesh2D
{
    /** The namespace for the functional module. **/
    namespace Functional
    {
    // namespace Function2D
    // {

    /*!
    * \addtogroup Functional
    * @{
    */

        template <unsigned rows, unsigned cols>
        using MatrixType = std::conditional_t<rows == 1 && cols == 1, double, Eigen::Matrix<double, rows, cols>>; ///< An alias for an Eigen::Matrix which reduces to a double when the matrix is of dimension 1x1.

        template <unsigned rows>
        using VectorType = MatrixType<rows, 1>; ///< An alias for a MatrixType with one column.

        /** 
         * The Function class is a wrapper class which stores a std::function describing the function value and another std::function describing the function derivative. \n\n
         * It describes in an abstract sense a function that maps elements from \f$\mathbb{R}^n \to \mathbb{R}^m\f$, for arbitrary postive integers \f$n, m\f$. \n\n
         * It is the base class for each of ScalarFunction1D, ScalarFunction2D, VectorFunction1D, VectorFunction2D, as well as Curve. It is the core object in the Functional module.
         * @tparam input_dim The dimension of the domain - must be a (strictly) positive integer type
         * @tparam output_dim The dimension of the range - must be a (strictly) positive integer type
        **/
        template <unsigned input_dim, unsigned output_dim>
        class Function
        {
        public:
            static_assert(input_dim > 0 && output_dim > 0, "ERROR - Function must have positive dimensions!\n");

            using InputType = VectorType<input_dim>; ///< The type of input the function expects.
            using OutputType = VectorType<output_dim>;  ///< The type of output the function value returns.
            using DerivativeType = MatrixType<output_dim, input_dim>; ///< The type of output the function derivative returns.

            using Value = std::function<OutputType(InputType)>; ///< An alias for a lambda function which defines the function value.
            using Derivative = std::function<DerivativeType(InputType)>; ///< An alias for a lambda function which defines the function derivative.

            //! A constructor for a Function with a value and a derivative.
            /*!
             * @param val The value of the function.
             * @param deriv The derivative of the function.
             */
            Function(const Value &val, const Derivative &deriv);

            //! A constructor for a Function which only stores a value.
            /*!
             * @param val The value of the function.
             */
            Function(const Value &val);

            //! A null constructor
            Function();

            //! A method to return the value of the function at a point x.
            virtual OutputType value(const InputType &x) const;

            //! A method to return the derivative of the function at a point x.
            virtual DerivativeType derivative(const InputType &x) const;

            //! A method to return the raw lambda describing the function value.
            const Value &get_value() const;

            //! A method to return the raw lambda describing the function derivative.
            const Derivative &get_derivative() const;

        private:
            Value __value; ///< A std::function describing the rule for the value of the function.
            Derivative __derivative; ///< A std::function describing the rule for the derivative of the function.

        // protected:
            //! The destructor is declared non-virtual and protected to prevent deletion of an inherited class through Function. This ensures the following code does not compile: VectorFunction2D *curve = new Curve(); delete curve; 
            // ~Function() {}
        };

        using ScalarFunction1D = Function<1, 1>; ///< An alias for a Function that maps R to R.
        using ScalarFunction2D = Function<2, 1>; ///< An alias for a Function that maps R^2 to R.
        using VectorFunction1D = Function<1, 2>; ///< An alias for a Function that maps R to R^2.
        using VectorFunction2D = Function<2, 2>; ///< An alias for a Function that maps R^2 to R^2.

        using ColVector = Eigen::Matrix<double, 2, 1>; ///< A column vector in two dimensions.
        using RowVector = Eigen::Matrix<double, 1, 2>; ///< A row vector in two dimensions.

        /** Function is a templated class describing an abstract function that maps elements from $R^n$ to $R^m$. It takes in two template parameters - the input dimension and output dimension. Function is a wrapper class which stores a std::function describing the function value and another std::function describing the function derivative.
        **/
        class Curve : public VectorFunction1D
        {
        public:
            Curve() {}
            Curve(const double t0, const double t1, const std::function<ColVector(double)> &val, const std::function<ColVector(double)> &deriv);

            double tmin;
            double tmax;

            ColVector value(const double t) const;
            ColVector derivative(const double t) const;

            ColVector tangent(const double t) const;
            ColVector normal(const double t) const;

            double derivative_norm(const double t) const;
            // double curvature(const double t) const;

        private:
            bool input_check(const double t) const;
            // std::function<double(double)> __curvature;
        };


        /** 
         * Given two functions \f$g:\mathbb{R}^n\to\mathbb{R}^k\f$ and \f$f:\mathbb{R}^k\to\mathbb{R}^m\f$, compose constructs a Function <n, m> 
         * for the composition \f$f\circ g:\mathbb{R}^n\to\mathbb{R}^m\f$.
         * @tparam input_dim (\f$n\f$ in description above) - The input dimension of the first function, \f$g\f$. Must be a strictly positive integer value.
         * @tparam intermediate_dim (\f$k\f$ in description above) - The output dimension of the first function, \f$g\f$, and input dimension of the second function, \f$f\f$. Must be a strictly positive integer value.
         * @tparam output_dim (\f$m\f$ in description above) - The output dimension of the second function, \f$f\f$. Must be a strictly positive integer value.
        **/
        template <unsigned input_dim, unsigned intermediate_dim, unsigned output_dim>
        const Function<input_dim, output_dim> compose(const Function<intermediate_dim, output_dim> &f, const Function<input_dim, intermediate_dim> &g);

        //@}

        // ---------------- Implementation of templated functions ----------------

        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::Function(const std::function<Function<input_dim, output_dim>::OutputType(Function<input_dim, output_dim>::InputType)> &val, const std::function<Function<input_dim, output_dim>::DerivativeType(Function<input_dim, output_dim>::InputType)> &deriv) : __value(val), __derivative(deriv) {}

        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::Function(const std::function<Function<input_dim, output_dim>::OutputType(Function<input_dim, output_dim>::InputType)> &val) : __value(val) {}

        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::Function() {}

        template <unsigned input_dim, unsigned output_dim>
        typename Function<input_dim, output_dim>::OutputType Function<input_dim, output_dim>::value(const Function<input_dim, output_dim>::InputType &x) const
        {
            return __value(x);
        }

        template <unsigned input_dim, unsigned output_dim>
        typename Function<input_dim, output_dim>::DerivativeType Function<input_dim, output_dim>::derivative(const Function<input_dim, output_dim>::InputType &x) const
        {
            return __derivative(x);
        }

        template <unsigned input_dim, unsigned output_dim>
        const typename Function<input_dim, output_dim>::Value &Function<input_dim, output_dim>::get_value() const
        {
            return __value;
        }

        template <unsigned input_dim, unsigned output_dim>
        const typename Function<input_dim, output_dim>::Derivative &Function<input_dim, output_dim>::get_derivative() const
        {
            return __derivative;
        }


        template <unsigned input_dim, unsigned intermediate_dim, unsigned output_dim>
        const Function<input_dim, output_dim> compose(const Function<intermediate_dim, output_dim> &f, const Function<input_dim, intermediate_dim> &g)
        {
            std::function<VectorType<output_dim>(VectorType<input_dim>)> val = [f, g](VectorType<input_dim> x) -> VectorType<output_dim>
            {
                return f.value(g.value(x));
            };
            std::function<MatrixType<output_dim, input_dim>(VectorType<input_dim>)> deriv = [f, g](VectorType<input_dim> x) -> MatrixType<output_dim, input_dim>
            {
                return f.derivative(g.value(x)) * g.derivative(x); // potential cast from Eigen::Matrix<double, 1, 1> to double
            };
            return Function<input_dim, output_dim>(val, deriv);
        }
    }
}

#endif
