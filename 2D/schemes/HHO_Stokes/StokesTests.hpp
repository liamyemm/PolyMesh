// Class to provide various test cases for the Stokes problem
//
//
// Author: Liam Yemm (liam.yemm@monash.edu)
//

#include "function.hpp"
#include "Mesh.hpp"
#include "math.hpp"
#include <Eigen/Dense>

#ifndef _STOKES_TESTS_HPP
#define _STOKES_TESTS_HPP

namespace PolyMesh2D
{
    class StokesSingularity
    {
    public:
        /// Initialise data
        StokesSingularity(double lambda, double alpha, double rot, Eigen::Vector2d x0);

        Functional::Function<2, 2> u() const;
        Functional::Function<2, 2> laplace_u() const;
        Functional::Function<2, 1> invcurl_u() const;
        Functional::Function<2, 1> p(double avg_pressure) const;

    private:
        double m_lambda;
        double m_alpha;
        double m_rot;
        Eigen::Vector2d m_x0;
    };

    // class SmoothStokesTestCircle
    // {
    // public:
    //     /// Initialise data
    //     SmoothStokesTestCircle() {}

    //     Functional::Function<2, 2> u() const
    //     {
    //         std::function<Eigen::Vector2d(Eigen::Vector2d)> u_lam = [](const Eigen::Vector2d x) -> Eigen::Vector2d
    //         {
    //             return (1.0 - x.squaredNorm()) * Eigen::Vector2d(-x(1), x(0));
    //         };
    //         std::function<Eigen::Matrix2d(Eigen::Vector2d)> Du_lam = [](const Eigen::Vector2d x) -> Eigen::Matrix2d
    //         {
    //             Eigen::Matrix2d mat = Eigen::Matrix2d::Zero();
    //             mat(0, 0) = 2.0 * x(0) * x(1);
    //             mat(0, 1) = -1.0 + x(0) * x(0) + 3.0 * x(1) * x(1);
    //             mat(1, 0) = 1.0 - 3.0 * x(0) * x(0) - x(1) * x(1);
    //             mat(1, 1) = -2.0 * x(0) * x(1);

    //             return mat;
    //         };
    //         return Functional::Function<2, 2>(u_lam, Du_lam);
    //     }
    //     Functional::Function<2, 2> laplace_u() const
    //     {
    //     }
    //     Functional::Function<2, 1> p() const
    //     {
    //         std::function<double(Eigen::Vector2d)> p_lam = [](const Eigen::Vector2d x) -> double
    //         {
    //             return std::exp(x(0) + x(1)) - std::pow(std::exp(1) - 1, 2);
    //         };
    //         std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Dp_lam = [](const Eigen::Vector2d x) -> PolyMesh2D::Functional::RowVector
    //         {
    //             return std::exp(x(0) + x(1)) * PolyMesh2D::Functional::RowVector(1, 1);
    //         };
    //         return Functional::Function<2, 1>(p_lam, Dp_lam);
    //     }
    // };

    class StokesTestSquare
    {
    public:
        /// Initialise data
        StokesTestSquare() {}

        Functional::Function<2, 2> u() const
        {
            std::function<Eigen::Vector2d(Eigen::Vector2d)> u_lam = [](const Eigen::Vector2d x) -> Eigen::Vector2d
            {
                return std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)) * Eigen::Vector2d(std::sin(Math::PI * x(0)) * std::cos(Math::PI * x(1)), -std::cos(Math::PI * x(0)) * std::sin(Math::PI * x(1)));
            };
            std::function<Eigen::Matrix2d(Eigen::Vector2d)> Du_lam = [](const Eigen::Vector2d x) -> Eigen::Matrix2d
            {
                Eigen::Matrix2d mat = Eigen::Matrix2d::Zero();
                mat(0, 0) = 2.0 * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)) * std::cos(Math::PI * x(0)) * std::cos(Math::PI * x(1));
                mat(0, 1) = std::pow(std::sin(Math::PI * x(0)) * std::cos(Math::PI * x(1)), 2) - std::pow(std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)), 2);
                mat(1, 0) = std::pow(std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)), 2) - std::pow(std::cos(Math::PI * x(0)) * std::sin(Math::PI * x(1)), 2);
                mat(1, 1) = -2.0 * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)) * std::cos(Math::PI * x(0)) * std::cos(Math::PI * x(1));

                return Math::PI * mat;
            };
            return Functional::Function<2, 2>(u_lam, Du_lam);
        }
        Functional::Function<2, 2> laplace_u() const
        {
            std::function<Eigen::Vector2d(Eigen::Vector2d)> laplace_u_lam = [](const Eigen::Vector2d x) -> Eigen::Vector2d
            {
                return -Math::PI * Math::PI * Eigen::Vector2d((1.0 - 2.0 * std::cos(2.0 * Math::PI * x(0))) * std::sin(2.0 * Math::PI * x(1)), -(1.0 - 2.0 * std::cos(2.0 * Math::PI * x(1))) * std::sin(2.0 * Math::PI * x(0)));
            };
            std::function<Eigen::Matrix2d(Eigen::Vector2d)> D_laplace_u_lam = [](const Eigen::Vector2d x) -> Eigen::Matrix2d
            {
                Eigen::Matrix2d mat = Eigen::Matrix2d::Zero();
                mat(0, 0) = 2.0 * std::sin(2.0 * Math::PI * x(0)) * std::sin(2.0 * Math::PI * x(1));
                mat(0, 1) = (1.0 - 2.0 * std::cos(2.0 * Math::PI * x(0))) * std::cos(2.0 * Math::PI * x(1));
                mat(1, 0) = -(1.0 - 2.0 * std::cos(2.0 * Math::PI * x(1))) * std::cos(2.0 * Math::PI * x(0));
                mat(1, 1) = -2.0 * std::sin(2.0 * Math::PI * x(0)) * std::sin(2.0 * Math::PI * x(1));

                return -2.0 * std::pow(Math::PI, 3) * mat;
            };
            return Functional::Function<2, 2>(laplace_u_lam, D_laplace_u_lam);
        }
        Functional::Function<2, 1> p() const
        {
            std::function<double(Eigen::Vector2d)> p_lam = [](const Eigen::Vector2d x) -> double
            {
                // return std::sin(2.0 * Math::PI * x(0)) * std::sin(2.0 * Math::PI * x(1));
                return 1E3 * (std::exp(x(0) + x(1)) - std::pow(std::exp(1) - 1, 2));
            };
            std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Dp_lam = [](const Eigen::Vector2d x) -> PolyMesh2D::Functional::RowVector
            {
                // return 2.0 * Math::PI * PolyMesh2D::Functional::RowVector(std::cos(2.0 * Math::PI * x(0)) * std::sin(2.0 * Math::PI * x(1)), std::sin(2.0 * Math::PI * x(0)) * std::cos(2.0 * Math::PI * x(1)));
                return 1E3 * (std::exp(x(0) + x(1)) * PolyMesh2D::Functional::RowVector(1, 1));
            };
            return Functional::Function<2, 1>(p_lam, Dp_lam);
        }
        Functional::Function<2, 2> grad_p() const
        {
            std::function<Eigen::Vector2d(Eigen::Vector2d)> grad_p_lam = [](const Eigen::Vector2d x) -> Eigen::Vector2d
            {
                return std::exp(x(0) + x(1)) * Eigen::Vector2d(1, 1);
            };
            std::function<Eigen::Matrix2d(Eigen::Vector2d)> D_grad_p_lam = [](const Eigen::Vector2d x) -> Eigen::Matrix2d
            {
                return std::exp(x(0) + x(1)) * Eigen::Matrix2d::Ones();
            };
            return Functional::Function<2, 2>(grad_p_lam, D_grad_p_lam);
        }


        Functional::Function<2, 2> source() const
        {
            std::function<Eigen::Vector2d(Eigen::Vector2d)> source_lam = [this](const Eigen::Vector2d x) -> Eigen::Vector2d
            {
                return -this->laplace_u().value(x) + this->p().derivative(x).transpose();
            };
            return Functional::Function<2, 2>(source_lam);
        }
    };
}

//@}

#endif
