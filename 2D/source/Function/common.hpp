#include "function.hpp"
#include "basis.hpp"
#include <array>
#include <Eigen/Dense>

#ifndef _COMMON_HPP
#define _COMMON_HPP

namespace PolyMesh2D
{
    namespace Functional
    {
        /*!
         * \addtogroup Functional
         * @{
         */

        /**
         * Constructs a Curve object describing an arclength parameterisation of a line segment starting at v0 and terminating at v1.
         **/
        const Curve StraightLine(const ColVector &v0, const ColVector &v1);

        /**
         * Constructs a scalar-valued monomial function in 2 dimensions.
         * @param powers is an array of two positive integer values storing the powers of each of the components of the monomial.
         **/
        const ScalarFunction2D Monomial(const std::array<size_t, 2> &powers);

        /**
         * Constructs a scalar-valued monomial function in 1 dimension.
         * @param power is a positive integer value which determines the power of the monomial.
         **/
        const ScalarFunction1D Monomial(const size_t power);

        /**
         * Constructs a vector-valued monomial function in 1 dimension.
         * Only one component of the vector is non zero. @param index determines which vector component is non-zero
         * @param power is a positive integer value which determines the power of the monomial.
         **/
        const VectorFunction1D VectorMonomial(const size_t power, const size_t index);

        /**
         * Constructs a vector-valued monomial function in 2 dimensions.
         * Only one component of the vector is non zero. @param index determines which vector component is non-zero
         * @param powers is an array of two positive integer values storing the powers of each of the components of the monomial.
         **/
        const VectorFunction2D VectorMonomial(const std::array<size_t, 2> &powers, const size_t index);

        /**
         * Generates a scalar-valued monomial basis in 2 dimensions.
         * @param transform describes a coordinate transformation under which the monomial is described. i.e. if the transformation is \f$\mathbf{y} = \mathbf{f}(\mathbf{x})\f$, then each basis function is a monomial in \f$\mathbf{y}\f$. transform is required to be affine (and non-degenerate) for the basis to span the same space as the polynomials in \f$\mathbf{x}\f$. No checks are performed to verify this.
         * @param degree determines the degree of the polynomial space.
         **/
        const ScalarBasis2D MonomialScalarBasis(const VectorFunction2D &transform, const size_t degree);

        /**
         * Generates a vector-valued monomial basis in 2 dimensions.
         * @param transform describes a coordinate transformation under which the monomial is described. i.e. if the transformation is \f$\mathbf{y} = \mathbf{f}(\mathbf{x})\f$, then each basis function is a monomial in \f$\mathbf{y}\f$. transform is required to be affine (and non-degenerate) for the basis to span the same space as the polynomials in \f$\mathbf{x}\f$. No checks are performed to verify this.
         * @param degree determines the degree of the polynomial space.
         **/
        const VectorBasis2D MonomialVectorBasis(const VectorFunction2D &transform, const size_t degree);

        /**
         * Generates a basis for the Raviart-Thomas-Nedelec space \f$\mathbb{RTN}^{k}(T) = \mathbb{P}^{k-1}(T)^d + \mathbf{x}\mathbb{P}^{k-1}(T)\f$ in 2 dimensions.
         * @param transform describes a coordinate transformation under which the basis is described. i.e. if the transformation is \f$\mathbf{y} = \mathbf{f}(\mathbf{x})\f$, then the basis for the the Raviart-Thomas-Nedelec is given by vector valued monomials in \f$\mathbf{y}\f$ plus \f$\mathbf{y}\f$ times scalar valued monomials in \f$\mathbf{y}\f$. transform is required to be affine (and non-degenerate) for the basis to span the same space as the Raviart-Thomas-Nedelec space in \f$\mathbf{x}\f$. No checks are performed to verify this.
         * @param degree determines the degree of the Raviart-Thomas-Nedelec space.
         **/
        const VectorBasis2D RaviartThomasNedelecBasis(const VectorFunction2D &transform, const size_t degree);

        /**
         * Generates a set of functions which span \f$\mathbb{P}^0(E) + \mathbb{P}^k(\mathbb{R}^2)^2\cdot\mathbf{n}_E \f$ on a given Curve \f$E\f$.
         * @param param The parameterisation of the Curve \f$E\f$.
         * @param degree The degree of the polynomial space.
         **/
        const ScalarBasis1D CurvedEdgeBasis(const Curve &param, const size_t degree);

        /**
         * Generates a set of functions which span \f$\mathbb{P}^0(E)^2 + \mathbb{P}^k(\mathbb{R}^2)^{2\times 2} \mathbf{n}_E \f$ on a given Curve \f$E\f$.
         * @param param The parameterisation of the Curve \f$E\f$.
         * @param degree The degree of the polynomial space.
         **/
        const VectorBasis1D CurvedEdgeVectorBasis(const Curve &param, const size_t degree);

        /**
         * Generates a scalar-valued shifted monomial basis in the parameter of a Curve. i.e. If the Curve is defined by \f$E=\{\gamma(t)\}\f$ then the basis is in the coordinate \f$t\f$.
         * @param param The parameterisation of the Curve \f$E\f$. Required to be affine (i.e. a straight line segment) in order for the space to be equivalent to the restriction of the polynomial space \f$\mathbb{P}^k(\mathbb{R}^2)\f$ to the Curve \f$E\f$.
         * @param degree The degree of the polynomial space.
         **/
        const ScalarBasis1D MonomialScalarBasis(const Curve &param, const size_t degree);

        /**
         * Generates a vector-valued shifted monomial basis in the parameter of a Curve. i.e. If the Curve is defined by \f$E=\{\gamma(t)\}\f$ then the basis is in the coordinate \f$t\f$.
         * @param param The parameterisation of the Curve \f$E\f$. Required to be affine (i.e. a straight line segment) in order for the space to be equivalent to the restriction of the polynomial space \f$\mathbb{P}^k(\mathbb{R}^2)\f$ to the Curve \f$E\f$.
         * @param degree The degree of the polynomial space.
         **/
        const VectorBasis1D MonomialVectorBasis(const Curve &param, const size_t degree);

        /**
         * Constructs a Curve which is the restriction to a sub-domain of another Curve.
         * @param curve The original Curve.
         * @param t0 The new tmin for the constructed Curve.
         * @param t1 The new tmax for the constructed Curve.
         **/
        const Curve restriction(const Curve &curve, const double t0, const double t1);

        /**
         * Uses a bisection method to find a point on a Curve which minimises the distance to a given point. This method assumes the signed curvature does not change sign so that
         * there exists at most one local maximum and one local minimum (given an infintite domain extension of the curve)
         * @param curve The Curve we are minimising the distance to.
         * @param point The point we are minimising the distance to.
         * @return The parameter value of the curve which minimises the distance to point.
         **/
        const double min_distance(const Curve &curve, const Eigen::Vector2d &point);

        /**
         * @brief Constructs a function object which returns the gradient of a 2D scalar function.
         *
         * @param func The function to compute the gradient of.
         * @return A Function object representing the gradient of the input function.
         */
        const Function<2, 2> gradient(const Function<2, 1> &func);

        /**
         * @brief Constructs a function object which returns the vector curl of a 2D scalar function.
         *
         * @param func The function to compute the curl of.
         * @return A Function object representing the curl of the input function.
         */
        const Function<2, 2> curl(const Function<2, 1> &func);

        /**
         * @brief Constructs a function object which returns the scalar curl of a 2D vector function.
         *
         * @param func The function to compute the curl of.
         * @return A Function object representing the curl of the input function.
         */
        const Function<2, 1> curl(const Function<2, 2> &func);

        /**
         * @brief Constructs a function object which returns the the dot product of a vector valued function with the normal of a given curve.
         *
         * @param func The function to compute dot product with.
         * @param param The curve whose normal we take the dot product with.
         * @return A Function object representing the dot product of a vector valued function with the normal.
         */
        const Function<1, 1> dot_n(const Function<1, 2> &func, const Curve &param);

        /**
         * @brief Constructs a function object which returns the the product of a vector valued function with the normal of a given curve.
         *
         * @param func The function to compute product with.
         * @param param The curve whose normal we take the product with.
         * @return A Function object representing the product of a vector valued function with the normal.
         */
        const Function<1, 2> times_n(const Function<1, 1> &func, const Curve &param);

        /**
         * @brief Computes the trace of a function along a given curve.
         *
         * This function calculates the trace of a given function along a specified curve. The trace represents the values of the function projected onto the curve.
         *
         * @tparam output_dim The number of output dimensions of the function.
         * @param func The function to compute the trace of.
         * @param param The curve along which the trace is computed.
         * @return A Function object representing the trace of the input function.
         */
        template <unsigned output_dim>
        const Function<1, output_dim> trace(const Function<2, output_dim> &func, const Curve &param)
        {
            std::function<VectorType<output_dim>(double)> val = [func, param](double x) -> VectorType<output_dim>
            {
                return func.value(param.value(x));
            };

            Function<1, output_dim> ret(val);

            auto old_poles = func.get_poles();
            std::vector<Pole<double>> new_poles;
            for (auto &pole : old_poles)
            {
                // we only consider poles occuring at the end points of the parameterisation
                if ((pole.location - param.value(param.tmin)).norm() < 1e-12)
                {
                    Pole<double> new_pole;
                    new_pole.location = param.tmin;
                    new_pole.order = pole.order;
                    new_poles.push_back(new_pole);
                }
                else if ((pole.location - param.value(param.tmax)).norm() < 1e-12)
                {
                    Pole<double> new_pole;
                    new_pole.location = param.tmax;
                    new_pole.order = pole.order;
                    new_poles.push_back(new_pole);
                }
                else
                {
                    // assert(false);
                }
            }

            ret.set_poles(new_poles);
            return ret;
        }

        /**
         * @brief Computes the Neumann trace of a function along a given curve.
         *
         * @tparam output_dim The number of output dimensions of the function.
         * @param func The function to compute the Neumann trace of.
         * @param param The curve along which the Neumann trace is computed.
         * @return A Function object representing the Neumann trace of the input function.
         */
        template <unsigned output_dim>
        const Function<1, output_dim> neumann_trace(const Function<2, output_dim> &func, const Curve &param, const Eigen::Vector2d &bias = Eigen::Vector2d::Zero())
        {
            std::function<VectorType<output_dim>(double)> val = [func, param, bias](double x) -> VectorType<output_dim>
            {
                return func.derivative(param.value(x) + bias) * param.normal(x);
            };

            Function<1, output_dim> ret(val);

            auto old_poles = func.get_poles();
            std::vector<Pole<double>> new_poles;
            for (auto &pole : old_poles)
            {
                // we only consider poles occuring at the end points of the parameterisation
                if ((pole.location - param.value(param.tmin)).norm() < 1e-12)
                {
                    Pole<double> new_pole;
                    new_pole.location = param.tmin;
                    new_pole.order = pole.order + 1;
                    new_poles.push_back(new_pole);
                }
                else if ((pole.location - param.value(param.tmax)).norm() < 1e-12)
                {
                    Pole<double> new_pole;
                    new_pole.location = param.tmax;
                    new_pole.order = pole.order + 1;
                    new_poles.push_back(new_pole);
                }
            }

            ret.set_poles(new_poles);
            return ret;
        }

        /**
         * @brief Constructs a function object which returns the scalar product of the gradients of two Function objects.
         *
         * @tparam output_dim The number of output dimensions of the input functions.
         * @param func1 The first function.
         * @param func2 The second function.
         * @return A Function object representing the scalar product of the gradients of two input functions.
         */
        template <unsigned output_dim>
        const Function<2, 1> gradient_product(const Function<2, output_dim> &func1, const Function<2, output_dim> &func2)
        {
            std::function<double(Eigen::Vector2d)> grad_1_dot_grad_2 = [func1, func2](const Eigen::Vector2d &x) -> double
            {
                return Math::scalar_product(func1.derivative(x), func2.derivative(x));
            };

            std::vector<Functional::Pole<Eigen::Vector2d>> poles1(func1.get_poles());
            std::vector<Functional::Pole<Eigen::Vector2d>> poles2(func2.get_poles());
            for (auto &pole : poles1)
                pole.order = pole.order + 1;
            for (auto &pole : poles2)
                pole.order = pole.order + 1;

            Function<2, 1> ret(grad_1_dot_grad_2);

            ret.set_poles(combine_poles(poles1, poles2));
            return ret;
        }

        /**
         * Extracts a component from a higher-dimensional function.
         *
         * This method takes a higher-dimensional function and extracts a single component from it, resulting in a new function with a reduced output dimension of 1. The component to extract is specified by the `index` parameter.
         *
         * @tparam input_dim The dimension of the input space.
         * @tparam output_dim The dimension of the output space of the original function.
         * @param func The original function from which to extract the component.
         * @param index The index of the component to extract.
         * @return A new function with a reduced output dimension of 1, representing the extracted component.
         *
         * @pre `index` must be less than `output_dim` of the original function.
         */
        template <unsigned input_dim, unsigned output_dim>
        const Function<input_dim, 1> component(const Function<input_dim, output_dim> &func, size_t index)
        {
            assert(index < output_dim);
            auto val = [func, index](const typename Function<input_dim, 1>::InputType &x) -> typename Function<input_dim, 1>::OutputType
            {
                return func.value(x)(index);
            };

            auto deriv = [func, index](const typename Function<input_dim, 1>::InputType &x) -> typename Function<input_dim, 1>::DerivativeType
            {
                return func.derivative(x).row(index);
            };

            Function<input_dim, 1> ret(val, deriv);
            ret.set_poles(func.get_poles());

            return ret;
        }

        //@}
    }
}
#endif
