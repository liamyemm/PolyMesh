#include "TestCase.hpp"
#include <iostream>

namespace PolyMesh2D
{
    namespace HHOPOISSON
    {
        TestCase::TestCase(unsigned id_test_case, char bdry_case) : m_id_test_case(id_test_case)
        {
            if (bdry_case == 'C')
            {
                std::function<Eigen::Vector2d(double)> circle_val = [](double t) -> Eigen::Vector2d
                { return Eigen::Vector2d(std::cos(t) - 0.25, std::sin(t) - 0.25); };
                std::function<Eigen::Vector2d(double)> circle_deriv = [](double t) -> Eigen::Vector2d
                { return Eigen::Vector2d(-std::sin(t), std::cos(t)); };

                bdry_param = Curve(0.0, 2.0 * Math::PI, circle_val, circle_deriv);

                std::function<double(Functional::ColVector)> LS = [](const Functional::ColVector &x) -> double
                {
                    return 1.0 - (std::pow(x(0) + 0.25, 2) + std::pow(x(1) + 0.25, 2));
                };

                std::function<Functional::RowVector(Functional::ColVector)> LS_grad = [](const Functional::ColVector &x) -> Functional::RowVector
                {
                    return -2.0 * Functional::RowVector(x(0), x(1));
                };

                level_set = ScalarFunction2D(LS, LS_grad);

                level_set_laplace = [](const Functional::ColVector &x) -> double
                {
                    return -4.0;
                };
            }
            else if(bdry_case == 'E')
            {
                double a = 0.8;

                std::function<Eigen::Vector2d(double)> ellipse_val = [a](double t) -> Eigen::Vector2d
                { return a * Eigen::Vector2d(std::cos(t) / std::sqrt(3) - std::sin(t), std::cos(t) / std::sqrt(3) + std::sin(t)); };
                std::function<Eigen::Vector2d(double)> ellipse_deriv = [a](double t) -> Eigen::Vector2d
                { return a * Eigen::Vector2d(-std::sin(t) / std::sqrt(3) - std::cos(t), -std::sin(t) / std::sqrt(3) + std::cos(t)); };

                bdry_param = Curve(0.0, 2.0 * Math::PI, ellipse_val, ellipse_deriv);

                std::function<double(Eigen::Vector2d)> LS = [a](const Eigen::Vector2d &x) -> double
                {
                    return a * a - (x(0) * x(0) + x(0) * x(1) + x(1) * x(1));
                };

                std::function<Functional::RowVector(Functional::ColVector)> LS_grad = [](const Functional::ColVector &x) -> Functional::RowVector
                {
                    return -Functional::RowVector(2.0 * x(0) + x(1), 2.0 * x(1) + x(0));
                };

                level_set = ScalarFunction2D(LS, LS_grad);

                level_set_laplace = [](const Functional::ColVector &x) -> double
                {
                    return -4.0;
                };
            }
            else
            {
                std::cout << "Error! Invalid boundary ID.\n";
                exit(1);
            }

            switch (m_id_test_case)
            {
            case 1:
                u = [](const double t) -> double
                {
                    return t;
                };
                Du = [](const double t) -> double
                {
                    return 1.0;
                };
                DDu = [](const double t) -> double
                {
                    return 0.0;
                };
                break;
            case 2:
                u = [](const double t) -> double
                {
                    return std::exp(t) - 1.0;
                };
                Du = [](const double t) -> double
                {
                    return std::exp(t);
                };
                DDu = [](const double t) -> double
                {
                    return std::exp(t);
                };
                break;
            case 3:
                u = [](const double t) -> double
                {
                    return std::sin(t);
                };
                Du = [](const double t) -> double
                {
                    return std::cos(t);
                };
                DDu = [](const double t) -> double
                {
                    return -std::sin(t);
                };
                break;
            default:
                break;
            }
        }

        Curve TestCase::get_boundary_param()
        {
            return bdry_param;
        }
        ScalarFunction2D TestCase::get_level_set()
        {
            return level_set;
        }

        // Solution
        ScalarFunction2D TestCase::sol()
        {
            Functional::ScalarFunction1D tmp(u, Du);
            return Functional::compose(tmp, level_set);
        }

        // source term
        ScalarFunction2D TestCase::src()
        {
            std::function<OutputType(InputType)> f = [&](const InputType &x) -> OutputType
            {
                return -level_set_laplace(x) * Du(level_set.value(x)) - level_set.derivative(x).squaredNorm() * DDu(level_set.value(x));
            };

            return ScalarFunction2D(f);
        }
    }
}