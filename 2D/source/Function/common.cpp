#include "common.hpp"
#include "math.hpp"

namespace PolyMesh2D
{
    namespace Functional
    {
        const Curve StraightLine(const ColVector &v0, const ColVector &v1)
        {
            double tmin = 0.0;
            double tmax = (v1 - v0).norm(); // arc length parameterisation
            std::function<ColVector(double)> val = [v0, v1](const double t) -> ColVector
            { return v0 + t * (v1 - v0) / (v1 - v0).norm(); };
            std::function<ColVector(double)> deriv = [v0, v1](const double t) -> ColVector
            { return (v1 - v0) / (v1 - v0).norm(); };

            return Curve(tmin, tmax, val, deriv);
        }

        const ScalarFunction2D Monomial(const std::array<size_t, 2> &powers)
        {
            std::function<double(ColVector)> val = [powers](const ColVector &x) -> double
            { return std::pow(x(0), powers[0]) * std::pow(x(1), powers[1]); };
            std::function<RowVector(ColVector)> deriv = [powers](const ColVector &x) -> RowVector
            {
                double x_deriv = (powers[0] == 0 ? 0.0 : powers[0] * std::pow(x(0), powers[0] - 1) * std::pow(x(1), powers[1]));
                double y_deriv = (powers[1] == 0 ? 0.0 : powers[1] * std::pow(x(0), powers[0]) * std::pow(x(1), powers[1] - 1));
                return RowVector(x_deriv, y_deriv);
            };
            return ScalarFunction2D(val, deriv);
        }

        const ScalarFunction1D Monomial(const size_t power)
        {
            std::function<double(double)> val = [power](const double x) -> double
            {
                return std::pow(x, power);
            };
            std::function<double(double)> deriv = [power](const double x) -> double
            {
                return (power == 0 ? 0.0 : power * std::pow(x, power - 1));
            };
            return ScalarFunction1D(val, deriv);
        }

        const VectorFunction2D VectorMonomial(const std::array<size_t, 2> &powers, const size_t d)
        {
            if (d > 1)
            {
                throw "ERROR - index out of bounds\n";
            }

            using InputType = VectorFunction2D::InputType;
            using OutputType = VectorFunction2D::OutputType;
            using DerivativeType = VectorFunction2D::DerivativeType;

            std::function<OutputType(InputType)> val = [powers, d](const InputType &x) -> OutputType
            {
                OutputType val = OutputType::Zero();
                val(d) = std::pow(x(0), powers[0]) * std::pow(x(1), powers[1]);
                return val;
            };
            std::function<DerivativeType(InputType)> deriv = [powers, d](const InputType &x) -> DerivativeType
            {
                double x_deriv = (powers[0] == 0 ? 0.0 : powers[0] * std::pow(x(0), powers[0] - 1) * std::pow(x(1), powers[1]));
                double y_deriv = (powers[1] == 0 ? 0.0 : powers[1] * std::pow(x(0), powers[0]) * std::pow(x(1), powers[1] - 1));

                Eigen::Matrix2d deriv = DerivativeType::Zero();
                deriv.row(d) = RowVector(x_deriv, y_deriv);
                return deriv;
            };
            return VectorFunction2D(val, deriv);
        }

        const ScalarBasis2D MonomialScalarBasis(const VectorFunction2D &transform, const size_t degree)
        {
            ScalarBasis2D basis;
            for (size_t l = 0; l <= degree; l++)
            {
                for (size_t i = 0; i <= l; i++)
                {
                    auto func = compose(Monomial({i, l - i}), transform);
                    basis.add_basis_function(func);
                } // for i
            }
            return basis;
        }

        const VectorBasis2D MonomialVectorBasis(const VectorFunction2D &transform, const size_t degree)
        {
            VectorBasis2D basis;
            for (size_t d = 0; d < 2; ++d)
            {
                for (size_t l = 0; l <= degree; l++)
                {
                    for (size_t i = 0; i <= l; i++)
                    {
                        auto func = compose(VectorMonomial({i, l - i}, d), transform);
                        basis.add_basis_function(func);
                    } // for i
                }
            }
            return basis;
        }

        const ScalarBasis1D CurvedEdgeBasis(const Curve &edge_param, const size_t degree)
        {
            Eigen::Vector2d xT = edge_param.value(0.5 * (edge_param.tmax + edge_param.tmin));
            double hT = edge_param.tmax - edge_param.tmin;

            std::function<VectorType<2>(VectorType<2>)> transform_val = [xT, hT](const VectorType<2> &x) -> VectorType<2>
            {
                return (1.0 / hT) * (x - xT);
            };

            std::function<MatrixType<2, 2>(VectorType<2>)> transform_deriv = [xT, hT](const VectorType<2> &x) -> MatrixType<2, 2>
            {
                return (1.0 / hT) * Eigen::Matrix2d::Identity();
            };

            VectorFunction2D transform(transform_val, transform_deriv);

            VectorBasis2D vector_basis = MonomialVectorBasis(transform, degree);
            ScalarBasis1D basis;
            std::function<double(double)> const_func_val = [](const double x) -> double
            { return 1.0; };
            std::function<double(double)> const_func_deriv = [](const double x) -> double
            { return 0.0; };

            ScalarFunction1D const_func(const_func_val, const_func_deriv);
            basis.add_basis_function(const_func);
            for (size_t i = 0; i < vector_basis.dimension(); ++i)
            {
                auto &vec_i = vector_basis.function(i).get_value();

                std::function<double(double)> func = [vec_i, edge_param](const double x) -> double
                {
                    return vec_i(edge_param.value(x)).dot(edge_param.normal(x));
                };

                basis.add_basis_function(func);
            }
            return basis;
        }

        const ScalarBasis1D MonomialScalarBasis(const Curve &edge_param, const size_t degree)
        {
            ScalarBasis1D basis;

            double tmin = edge_param.tmin;
            double tmax = edge_param.tmax;

            std::function<double(double)> transform_val = [tmin, tmax](const double t) -> double
            {
                return (t - tmin) / (tmax - tmin);
            };

            std::function<double(double)> transform_deriv = [tmin, tmax](const double t) -> double
            {
                return 1.0 / (tmax - tmin);
            };

            ScalarFunction1D transform(transform_val, transform_deriv);

            for (size_t i = 0; i <= degree; ++i)
            {
                // basis.add_basis_function(Monomial(i));
                basis.add_basis_function(compose(Monomial(i), transform));
            }
            return basis;
        }

        const VectorBasis2D RaviartThomasNedelecBasis(const VectorFunction2D &transform, const size_t degree)
        {
            VectorBasis2D basis;
            for (size_t d = 0; d < 2; ++d)
            {
                for (size_t l = 0; l <= degree - 1; l++)
                {
                    for (size_t i = 0; i <= l; i++)
                    {
                        auto func = compose(VectorMonomial({i, l - i}, d), transform);
                        basis.add_basis_function(func);
                    } // for i
                }
            }
            for (size_t i = 0; i <= degree - 1; i++)
            {
                auto func = compose(Monomial({i, degree - 1 - i}), transform);
                basis.add_basis_function(func * transform);
            } 
            return basis;
        }

        double bisection_method(const std::function<double(double)> &f, const double t0, const double t1)
        {
            assert(f(t0) * f(t1) < 0.0);
            assert(t0 < t1);

            double tolerance = 1E-15;
            int max_its = 1000;

            double a = t0;
            double b = t1;

            int its = 0;
            while (its < max_its)
            {
                double c = 0.5 * (a + b);
                if (f(c) == 0.0 || 0.5 * (b - a) < tolerance)
                {
                    return c;
                }
                if(Math::sgn(f(c)) == Math::sgn(f(a)))
                {
                    a = c;
                }
                else
                {
                    b = c;
                }
                ++its;
            }
            assert(false && "Error! bisection_method failed to converge.\n");
            return 0.0;
        }

        const double min_distance(const Curve &curve, const Eigen::Vector2d &point)
        {
            assert((curve.value(curve.tmax) - curve.value(curve.tmin)).norm() < 1E-15); // must be a closed curve

            std::vector<double> t_vals;

            std::function<double(double)> find_roots_of = [&curve, &point](const double t) -> double
            {
                return (point - curve.value(t)).dot(curve.derivative(t));
            };

            unsigned n_partitions = 10; // probably overkill, but I am not certain.
            double delta = (curve.tmax - curve.tmin) / n_partitions;
            for(unsigned i = 0; i < n_partitions; ++i)
            {
                double t_start = curve.tmin + i * delta;
                double t_end = curve.tmin + (i + 1) * delta;

                if(std::abs(find_roots_of(t_start)) < 1E-15)
                {
                    t_vals.push_back(t_start);
                }
                else if (Math::sgn(find_roots_of(t_start)) * Math::sgn(find_roots_of(t_end)) < 0)
                {
                    t_vals.push_back(bisection_method(find_roots_of, t_start, t_end));
                }
            }

            assert(t_vals.size() >= 2); // should have at least two local optima
            // assert(t_vals.size() <= 4); // should have at most four local optima

            double t = t_vals[0];

            for(size_t i = 1; i < t_vals.size(); ++i)
            {
                if((curve.value(t_vals[i]) - point).norm() < (curve.value(t) - point).norm())
                {
                    t = t_vals[i]; // find the minimiser
                }
            }

            return t;
        }


        const Curve restriction(const Curve &curve, const double t0, const double t1)
        {
            double tmin = std::min(t0, t1);
            double tmax = std::max(t0, t1);

            if((curve.value(curve.tmax) - curve.value(curve.tmin)).norm() < 1E-15) // closed curve
            {
                bool tmax_fourth_quad = (tmax > curve.tmin + 0.75 * (curve.tmax - curve.tmin));
                bool tmin_first_quad = (tmin < curve.tmin + 0.25 * (curve.tmax - curve.tmin));
                if(tmin_first_quad && tmax_fourth_quad)
                {
                    tmin += curve.tmax; // assumes periodic. will cause issues otherwise
                    double tmp = tmax;
                    tmax = tmin;
                    tmin = tmp;
                }
            }

            return Curve(tmin, tmax, curve.get_value(), curve.get_derivative());
        }
    }
}