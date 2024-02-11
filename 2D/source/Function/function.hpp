// Provides wrapper classes and methods for abstract functions that map R^n to R^m
//
// Author: Liam Yemm (liam.yemm@monash.edu)

#ifndef _FUNCTION_HPP
#define _FUNCTION_HPP

#include <functional> // std::function
#include "math.hpp"
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
        /*!
         * \addtogroup Functional
         * @{
         */

        template <unsigned rows, unsigned cols>
        using MatrixType = std::conditional_t<rows == 1 && cols == 1, double, Eigen::Matrix<double, rows, cols>>; ///< An alias for an Eigen::Matrix which reduces to a double when the matrix is of dimension 1x1.

        template <unsigned rows>
        using VectorType = MatrixType<rows, 1>; ///< An alias for a MatrixType with one column.

        template <typename PointType>
        struct Pole
        {
            Pole() {}
            Pole(const PointType &loc, const double &ord) : location(loc), order(ord) {}
            PointType location;
            double order;
        };

        template <typename PointType>
        std::vector<Pole<PointType>> combine_poles(const std::vector<Pole<PointType>> &poles1, const std::vector<Pole<PointType>> &poles2)
        {
            std::vector<Pole<PointType>> combinedPoles;

            // Copy poles from the first vector
            combinedPoles.insert(combinedPoles.end(), poles1.begin(), poles1.end());

            // Check poles from the second vector
            for (const auto &pole2 : poles2)
            {
                bool found = false;

                // Search for a pole with the same location
                for (auto &pole1 : combinedPoles)
                {
                    PointType pole_diff(pole1.location - pole2.location);
                    if (Math::norm(pole_diff) < 1E-12)
                    {
                        // Pole with the same location found, combine the orders
                        pole1.order = pole1.order + pole2.order;
                        found = true;
                        break;
                    }
                }

                // If no pole with the same location was found, add it to the combined poles
                if (!found)
                {
                    combinedPoles.push_back(pole2);
                }
            }

            return combinedPoles;
        }

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

            using InputType = VectorType<input_dim>;                  ///< The type of input the function expects.
            using OutputType = VectorType<output_dim>;                ///< The type of output the function value returns.
            using DerivativeType = MatrixType<output_dim, input_dim>; ///< The type of output the function derivative returns.

            using Value = std::function<OutputType(InputType)>;          ///< An alias for a std::function which maps InputType to OutputType.
            using Derivative = std::function<DerivativeType(InputType)>; ///< An alias for a std::function which maps InputType to DerivativeType.

            virtual ~Function();

            //! A constructor for a Function with a value and a derivative.
            /*!
             * @param val The value of the function.
             * @param deriv The derivative of the function.
             */
            explicit Function(const Value &val, const Derivative &deriv);

            //! A constructor for a Function which only stores a value.
            /*!
             * @param val The value of the function.
             */
            explicit Function(const Value &val);

            //! Default constructor
            explicit Function();

            //! A method to return the value of the function at a point x.
            virtual OutputType value(const InputType &x) const;

            //! A method to return the derivative of the function at a point x.
            virtual DerivativeType derivative(const InputType &x) const;

            /**
             * Set the poles of the function.
             *
             * @param poles The vector of poles to be set.
             */
            void set_poles(const std::vector<Pole<InputType>> &poles);

            /**
             * Get the poles of the function.
             *
             * @return The vector of poles.
             */
            std::vector<Pole<InputType>> get_poles() const;

            /**
             * Get a specific pole of the function by index.
             *
             * @param index The index of the pole to retrieve.
             * @return The pole at the specified index.
             */
            Pole<InputType> pole(size_t index) const;

            /**
             * Add a pole to the function.
             *
             * @param pole The pole to be added.
             */
            void add_pole(const Pole<InputType> &pole);

            /**
             * Check if the function is differentiable.
             *
             * @return True if the function is differentiable, false otherwise.
             */
            bool differentiable() const;

        private:
            Value __value;           ///< A std::function describing the rule for the value of the function.
            Derivative __derivative; ///< A std::function describing the rule for the derivative of the function.

            bool __differentiable;
            std::vector<Pole<InputType>> poles;

        };

        /**
         * Overloaded operator to add two functions.
         * @tparam input_dim The number of input dimensions.
         * @tparam output_dim The number of output dimensions.
         * @param f1 The first function to be added.
         * @param f2 The second function to be added.
         * @return A new Function object which is the sum of f1 and f2.
         **/
        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim> operator+(const Function<input_dim, output_dim> &f1, const Function<input_dim, output_dim> &f2)
        {
            auto sum_value = [f1, f2](const typename Function<input_dim, output_dim>::InputType &x) -> typename Function<input_dim, output_dim>::OutputType
            {
                return f1.value(x) + f2.value(x);
            };

            if (f1.differentiable() && f2.differentiable()) // Check if both f1 and f2 have available derivatives
            {
                auto sum_derivative = [f1, f2](const typename Function<input_dim, output_dim>::InputType &x) -> typename Function<input_dim, output_dim>::DerivativeType
                {
                    return f1.derivative(x) + f2.derivative(x);
                };

                return Function<input_dim, output_dim>(sum_value, sum_derivative);
            }
            else
            {
                return Function<input_dim, output_dim>(sum_value);
            }

            // #warning "The + operator does not copy the poles of the input functions"
        }

        // /**
        //  * Overloaded operator to subtract two functions.
        //  * @tparam input_dim The number of input dimensions.
        //  * @tparam output_dim The number of output dimensions.
        //  * @param f1 The first function to be subtracted.
        //  * @param f2 The second function to be subtracted.
        //  * @return A new Function object which is the f1 - f2.
        //  **/
        // template <unsigned input_dim, unsigned output_dim>
        // Function<input_dim, output_dim> operator-(const Function<input_dim, output_dim> &f1, const Function<input_dim, output_dim> &f2)
        // {
        //     auto subtract_value = [f1, f2](const typename Function<input_dim, output_dim>::InputType &x) -> typename Function<input_dim, output_dim>::OutputType
        //     {
        //         return f1.value(x) - f2.value(x);
        //     };

        //     std::vector<Pole<VectorType<input_dim>>> poles = f1.poles;
        //     poles.insert(poles.end(), f2.poles.begin(), f2.poles.end());

        //     for (size_t i = 0; i < poles.size(); ++i)
        //     {
        //         for (size_t j = i + 1; j < poles.size(); ++j)
        //         {
        //             VectorType<input_dim> pole_diff(poles[i].location - poles[j].location);
        //             if (Math::norm(pole_diff) < 1E-15)
        //             {
        //                 poles[i].order = std::max(poles[i].order, poles[j].order);
        //                 poles.erase(poles.begin() + j);
        //                 if (j == i + 1)
        //                 {
        //                     --i;
        //                 }
        //                 break;
        //             }
        //         }
        //     }

        //     if (f1.differentiable() && f2.differentiable()) // Check if both f1 and f2 have available derivatives
        //     {
        //         auto subtract_derivative = [f1, f2](const typename Function<input_dim, output_dim>::InputType &x) -> typename Function<input_dim, output_dim>::DerivativeType
        //         {
        //             return f1.derivative(x) - f2.derivative(x);
        //         };

        //         Function<input_dim, output_dim> ret(subtract_value, subtract_derivative);
        //         ret.poles = poles;
        //         return ret;
        //     }
        //     else
        //     {
        //         Function<input_dim, output_dim> ret(subtract_value);
        //         ret.poles = poles;
        //         return ret;
        //     }
        // }

        /**
         * Overloaded operator to multiply a function by a scalar valued function.
         * @tparam input_dim The number of input dimensions.
         * @tparam output_dim The number of output dimensions.
         * @param f1 The first function to be multiplied - must be scalar valued.
         * @param f2 The second function to be multiplied.
         * @return A new Function object which is the product of f1 and f2.
         **/
        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim> operator*(const Function<input_dim, 1> &f1, const Function<input_dim, output_dim> &f2)
        {
            auto product_value = [f1, f2](const typename Function<input_dim, output_dim>::InputType &x) -> typename Function<input_dim, output_dim>::OutputType
            {
                return f1.value(x) * f2.value(x);
            };

            // std::vector<Pole<VectorType<input_dim>>> poles = f1.poles;
            // poles.insert(poles.end(), f2.poles.begin(), f2.poles.end());

            // for (size_t i = 0; i < poles.size(); ++i)
            // {
            //     for (size_t j = i + 1; j < poles.size(); ++j)
            //     {
            //         VectorType<input_dim> pole_diff(poles[i].location - poles[j].location);
            //         if (Math::norm(pole_diff) < 1E-15)
            //         {
            //             poles[i].order += poles[j].order;
            //             poles.erase(poles.begin() + j);
            //             if (j == i + 1)
            //             {
            //                 --i;
            //             }
            //             break;
            //         }
            //     }
            // }

            if (f1.differentiable() && f2.differentiable()) // Check if both f1 and f2 have available derivatives
            {
                auto product_derivative = [f1, f2](const typename Function<input_dim, output_dim>::InputType &x) -> typename Function<input_dim, output_dim>::DerivativeType
                {
                    return f1.value(x) * f2.derivative(x) + f2.value(x) * f1.derivative(x);
                };

                Function<input_dim, output_dim> ret(product_value, product_derivative);
                ret.set_poles(combine_poles(f1.get_poles(), f2.get_poles()));
                return ret;
            }
            else
            {
                Function<input_dim, output_dim> ret(product_value);
                ret.set_poles(combine_poles(f1.get_poles(), f2.get_poles()));
                return ret;
            }
        }

        // // included so that scalar_product compiles with vector valued functions. Should never actually be called
        // template <unsigned input_dim>
        // Function<input_dim, 1> operator*(const Function<input_dim, 2> &f1, const Function<input_dim, 2> &f2)
        // {
        //     throw std::runtime_error("Operator * cannot be called with two vector-valued functions");
        //     return Function<input_dim, 1>();
        // }

        using ScalarFunction1D = Function<1, 1>; ///< An alias for a Function that maps R to R.
        using ScalarFunction2D = Function<2, 1>; ///< An alias for a Function that maps R^2 to R.
        using VectorFunction1D = Function<1, 2>; ///< An alias for a Function that maps R to R^2.
        using VectorFunction2D = Function<2, 2>; ///< An alias for a Function that maps R^2 to R^2.

        using ColVector = Eigen::Matrix<double, 2, 1>; ///< A column vector in two dimensions.
        using RowVector = Eigen::Matrix<double, 1, 2>; ///< A row vector in two dimensions.

        /** A Curve is a function which maps values from \f$[t_0,t_1]\f$ to \f$R^2\f$. It is a sub-class of VectorFunction1D and stores the extra functions tangent, normal, and derivative_norm. Curve also stores parameters tmin and tmax which determine the end points of the curve.
         **/
        class Curve : public VectorFunction1D
        {
        public:
            //! A null constructor for a Curve.
            Curve() {}

            //! The Curve constructor.
            /*!
             * @param tmin The start point of the Curve in parameter space.
             * @param tmax The end point of the Curve in parameter space.
             * @param val The value of the function.
             * @param deriv The derivative of the function.
             */
            Curve(const double t0, const double t1, const std::function<ColVector(double)> &val, const std::function<ColVector(double)> &deriv);

            //! The start point of the Curve.
            double tmin;

            //! The end point of the Curve.
            double tmax;

            //! A method to return the value of the Curve at the parameter t.
            ColVector value(const double &t) const override;

            //! A method to return the derivative of the Curve at the parameter t.
            ColVector derivative(const double &t) const override;

            //! A method to return the tangent vector of the Curve at the parameter t.
            ColVector tangent(const double t) const;

            //! A method to return the normal vector of the Curve at the parameter t.
            ColVector normal(const double t) const;

            //! A method to return the norm of the derivative of the Curve at the parameter t.
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

        /**
         * @brief Computes the scalar product of two functions.
         *
         * This function computes the scalar product of two scalar-valued functions `f`
         * and `g`, which must have the same input dimension `input_dim`. The result is a scalar-valued function with input
         * dimension `input_dim`.
         *
         * @tparam input_dim The input dimension of the functions.
         * @param f The first function.
         * @param g The second function.
         * @return A scalar-valued function representing the scalar product of `f` and `g`.
         */
        template <unsigned input_dim>
        const Function<input_dim, 1> scalar_product(const Function<input_dim, 1> &f, const Function<input_dim, 1> &g);

        /**
         * @brief Computes the scalar product of two functions.
         *
         * This function computes the scalar product of two vector-valued functions `f`
         * and `g`, which must have the same input dimension `input_dim`. The result is a scalar-valued function with input
         * dimension `input_dim`.
         *
         * @tparam input_dim The input dimension of the functions.
         * @param f The first function.
         * @param g The second function.
         * @return A scalar-valued function representing the scalar product of `f` and `g`.
         */
        template <unsigned input_dim>
        const Function<input_dim, 1> scalar_product(const Function<input_dim, 2> &f, const Function<input_dim, 2> &g);

        //@}

        // ---------------- Implementation of templated functions ----------------
        
        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::~Function()
        {
            // Destructor implementation
        }

        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::Function(const std::function<Function<input_dim, output_dim>::OutputType(Function<input_dim, output_dim>::InputType)> &val, const std::function<Function<input_dim, output_dim>::DerivativeType(Function<input_dim, output_dim>::InputType)> &deriv) : __value(val), __derivative(deriv), __differentiable(true) {}

        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::Function(const std::function<Function<input_dim, output_dim>::OutputType(Function<input_dim, output_dim>::InputType)> &val) : __value(val), __differentiable(false) {}

        template <unsigned input_dim, unsigned output_dim>
        Function<input_dim, output_dim>::Function() : __differentiable(false) {}

        template <unsigned input_dim, unsigned output_dim>
        typename Function<input_dim, output_dim>::OutputType Function<input_dim, output_dim>::value(const Function<input_dim, output_dim>::InputType &x) const
        {
            return __value(x);
        }

        template <unsigned input_dim, unsigned output_dim>
        typename Function<input_dim, output_dim>::DerivativeType Function<input_dim, output_dim>::derivative(const Function<input_dim, output_dim>::InputType &x) const
        {
            if (!differentiable())
            {
                // Throw an error indicating that the derivative is not available
                throw std::runtime_error("Call to derivative not available");
            }
            return __derivative(x);
        }

        template <unsigned input_dim, unsigned output_dim>
        void Function<input_dim, output_dim>::set_poles(const std::vector<Pole<InputType>> &poles)
        {
            this->poles = poles;
        }

        template <unsigned input_dim, unsigned output_dim>
        std::vector<Pole<typename Function<input_dim, output_dim>::InputType>> Function<input_dim, output_dim>::get_poles() const
        {
            return poles;
        }

        template <unsigned input_dim, unsigned output_dim>
        Pole<typename Function<input_dim, output_dim>::InputType> Function<input_dim, output_dim>::pole(size_t index) const
        {
            return poles.at(index);
        }

        template <unsigned input_dim, unsigned output_dim>
        void Function<input_dim, output_dim>::add_pole(const Pole<InputType> &pole)
        {
            poles.push_back(pole);
        }

        template <unsigned input_dim, unsigned output_dim>
        bool Function<input_dim, output_dim>::differentiable() const
        {
            return __differentiable;
        }

        template <unsigned input_dim, unsigned intermediate_dim, unsigned output_dim>
        const Function<input_dim, output_dim> compose(const Function<intermediate_dim, output_dim> &f, const Function<input_dim, intermediate_dim> &g)
        {
            std::function<VectorType<output_dim>(VectorType<input_dim>)> val = [f, g](VectorType<input_dim> x) -> VectorType<output_dim>
            {
                return f.value(g.value(x));
            };

            if (f.differentiable() && g.differentiable()) // Check if both f and g have available derivatives
            {
                std::function<MatrixType<output_dim, input_dim>(VectorType<input_dim>)> deriv = [f, g](VectorType<input_dim> x) -> MatrixType<output_dim, input_dim>
                {
                    return f.derivative(g.value(x)) * g.derivative(x); // potential cast from Eigen::Matrix<double, 1, 1> to double
                };

                Function<input_dim, output_dim> ret(val, deriv);
                return ret;
            }
            else
            {
                Function<input_dim, output_dim> ret(val);
                return ret;
            }
        }

        template <unsigned input_dim>
        const Function<input_dim, 1> scalar_product(const Function<input_dim, 1> &f, const Function<input_dim, 1> &g)
        {
            return f * g;
        }

        template <unsigned input_dim>
        const Function<input_dim, 1> scalar_product(const Function<input_dim, 2> &f, const Function<input_dim, 2> &g)
        {
            std::function<double(VectorType<input_dim>)> val = [f, g](VectorType<input_dim> x) -> double
            {
                return Math::scalar_product(f.value(x), g.value(x));
            };

            // if (f.differentiable() && g.differentiable()) // Check if both f and g have available derivatives
            // {
            //     std::function<MatrixType<1, input_dim>(VectorType<input_dim>)> deriv = [f, g](VectorType<input_dim> x) -> MatrixType<1, input_dim>
            //     {
            //         return f.derivative(x).trace() * g.value(x).transpose() + g.derivative(x).trace() * f.value(x).transpose();
            //     };

            //     Function<input_dim, 1> ret(val, deriv);
            //     ret.set_poles(combine_poles(f.get_poles(), g.get_poles()));
            //     return ret;
            // }
            // else
            // {
                Function<input_dim, 1> ret(val);
                ret.set_poles(combine_poles(f.get_poles(), g.get_poles()));
                return ret;
            // }
        }
    }
}

#endif
