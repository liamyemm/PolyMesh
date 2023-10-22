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
                mat(0, 1) = 1.0 / Math::PI + std::pow(std::sin(Math::PI * x(0)) * std::cos(Math::PI * x(1)), 2) - std::pow(std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)), 2);
                mat(1, 0) = 1.0 / Math::PI + std::pow(std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)), 2) - std::pow(std::cos(Math::PI * x(0)) * std::sin(Math::PI * x(1)), 2);
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
        Functional::Function<2, 1> p(double pressure_const) const
        {
            std::function<double(Eigen::Vector2d)> p_lam = [pressure_const](const Eigen::Vector2d x) -> double
            {
                // return std::sin(2.0 * Math::PI * x(0)) * std::sin(2.0 * Math::PI * x(1));
                return (x(0) - 0.5) * (x(1) - 0.5) * (x(1) - 0.5)  - pressure_const;
            };
            std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Dp_lam = [](const Eigen::Vector2d x) -> PolyMesh2D::Functional::RowVector
            {
                // return 2.0 * Math::PI * PolyMesh2D::Functional::RowVector(std::cos(2.0 * Math::PI * x(0)) * std::sin(2.0 * Math::PI * x(1)), std::sin(2.0 * Math::PI * x(0)) * std::cos(2.0 * Math::PI * x(1)));
                return PolyMesh2D::Functional::RowVector((x(1) - 0.5) * (x(1) - 0.5), 2.0 * (x(0) - 0.5) * (x(1) - 0.5));
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
                return -this->laplace_u().value(x) + this->p(0.0).derivative(x).transpose();
            };
            return Functional::Function<2, 2>(source_lam);
        }
    };

    class Circle
    {
    public:
        Circle(Eigen::Vector2d O, double R) : m_O(O), m_R(R) {}

        Functional::Curve param() const
        {
            Eigen::Vector2d O = m_O;
            double R = m_R;
            std::function<Eigen::Vector2d(double)> circle_val = [O, R](double t) -> Eigen::Vector2d
            { return R * Eigen::Vector2d(std::cos(t), std::sin(t)) + O; };
            std::function<Eigen::Vector2d(double)> circle_deriv = [R](double t) -> Eigen::Vector2d
            { return R * Eigen::Vector2d(-std::sin(t), std::cos(t)); };

            return Functional::Curve(0.0, 2.0 * Math::PI, circle_val, circle_deriv);
        }

        Functional::Function<2, 1> level_set() const
        {
            Eigen::Vector2d O = m_O;
            double R = m_R;
            std::function<double(PolyMesh2D::Functional::ColVector)> LS = [O, R](const PolyMesh2D::Functional::ColVector &x) -> double
            {
                return R * R - (x - O).squaredNorm();
            };

            std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> LS_grad = [O](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
            {
                return -2.0 * (x - O).transpose();
            };
            return Functional::Function<2, 1>(LS, LS_grad);
        }

        Functional::Function<2, 2> u() const
        {
            Eigen::Vector2d O = m_O;
            double R = m_R;
            std::function<Eigen::Vector2d(Eigen::Vector2d)> u_lam = [O, R](const Eigen::Vector2d point) -> Eigen::Vector2d
            {
                Eigen::Vector2d y = point - O;
                double r = y.norm();

                Eigen::Vector2d vec1(y(0) * y(0) - y(1) * y(1), 2 * y(0) * y(1));
                Eigen::Vector2d vec2(std::log(r / R), 0);

                return 0.5 * ((R * R - r * r) / std::pow(r, 4)) * vec1 + vec2;
            };
            std::function<Eigen::Matrix2d(Eigen::Vector2d)> Du_lam = [O, R](const Eigen::Vector2d point) -> Eigen::Matrix2d
            {
                Eigen::Vector2d y = point - O;
                double rsq = y.squaredNorm();

                Eigen::Matrix2d mat1 = Eigen::Matrix2d::Zero();
                mat1(0, 0) = y(0);
                mat1(0, 1) = y(1);
                mat1(1, 0) = y(1);
                mat1(1, 1) = -y(0);

                Eigen::Matrix2d mat2 = Eigen::Matrix2d::Zero();
                mat2(0, 0) = -y(1);
                mat2(0, 1) = y(0);
                mat2(1, 0) = y(0);
                mat2(1, 1) = y(1);

                Eigen::Matrix2d mat3 = Eigen::Matrix2d::Zero();
                mat3(0, 1) = 1.0;

                return ((rsq - R * R) * (y(0) * y(0) - y(1) * y(1)) / std::pow(rsq, 3)) * mat1 - (2.0 * y(0) * y(1) * R * R / std::pow(rsq, 3)) * mat2 + (2.0 * y(1) / rsq) * mat3;
            };
            return Functional::Function<2, 2>(u_lam, Du_lam);
        }

        Functional::Function<2, 2> laplace_u() const
        {
            Eigen::Vector2d O = m_O;
            std::function<Eigen::Vector2d(Eigen::Vector2d)> lap_u_lam = [O](const Eigen::Vector2d point) -> Eigen::Vector2d
            {
                Eigen::Vector2d y = point - O;
                double rsq = y.squaredNorm();

                return (2.0 / (rsq * rsq)) * Eigen::Vector2d(y(0) * y(0) - y(1) * y(1), 2 * y(0) * y(1));
            };
            std::function<Eigen::Matrix2d(Eigen::Vector2d)> Dlap_u_lam = [O](const Eigen::Vector2d point) -> Eigen::Matrix2d
            {
                Eigen::Vector2d y = point - O;
                double rsq = y.squaredNorm();

                double on_diag = y(0) * (y(0) * y(0) - 3.0 * y(1) * y(1));
                double off_diag = y(1) * (y(1) * y(1) - 3.0 * y(0) * y(0));

                Eigen::Matrix2d mat = Eigen::Matrix2d::Zero();
                mat(0, 0) = -on_diag;
                mat(0, 1) = off_diag;
                mat(1, 0) = off_diag;
                mat(1, 1) = on_diag;

                return (4.0 / (rsq * rsq * rsq)) * mat;
            };
            return Functional::Function<2, 2>(lap_u_lam, Dlap_u_lam);
        }

        Functional::Function<2, 1> p(double pressure_const) const
        {
            Eigen::Vector2d O = m_O;

            std::function<double(Eigen::Vector2d)> value = [O, pressure_const](const Eigen::Vector2d &point) -> double
            {
                Eigen::Vector2d y = point - O;
                return -pressure_const - 2.0 * y(0) / y.squaredNorm();
            };

            std::function<Eigen::RowVector2d(Eigen::Vector2d)> deriv = [O](const Eigen::Vector2d &point) -> Eigen::RowVector2d
            {
                Eigen::Vector2d y = point - O;
                double rsq = y.squaredNorm();

                return (2.0 / (rsq * rsq)) * Eigen::RowVector2d(y(0) * y(0) - y(1) * y(1), 2 * y(0) * y(1));
            };

            return Functional::Function<2, 1>(value, deriv);
        }

        bool enrich_cell(CurvedMesh::Cell *c, double cut_off)
        {
            if(cut_off < 0.0)
            {
                return true;
            }
            if((c->center_mass() - m_O).norm() < cut_off + m_R) // within distance cut_off from circle edge
            {
                return true;
            }
            return false;
        }

        bool enrich_edge(CurvedMesh::Edge *e, double cut_off)
        {
            if(cut_off < 0.0)
            {
                return true;
            }
            for(auto & c : e->get_cells())
            {
                if(this->enrich_cell(c, cut_off))
                {
                    return true;
                }
            }
            return false;
        }

        bool p_linearly_independent(CurvedMesh::Edge *e) const
        {
            if (std::abs(e->vertex(0)->coords()(0) - e->vertex(1)->coords()(0)) < 1E-14)
            {
                // x const on edge
                if (std::abs(e->vertex(0)->coords()(0) - m_O(0)) < 1E-14)
                {
                    return false;
                }
            }
            if (!(e->is_straight()))
            {
                // edge is curved, test if it is an edge of this circle

                // if (std::abs((e->vertex(0)->coords() - m_O).norm() - m_R) < 1E-12)
                // {
                //     assert(std::abs((e->vertex(1)->coords() - m_O).norm() - m_R) < 1E-12);
                //     return false;
                // }
                return false;
            }

            return true;
        }

        bool u_linearly_independent(CurvedMesh::Edge *e) const
        {
            // should only be called if edge_degree >= 1
            if (!(e->is_straight()))
            {
                // edge is on a circle, test if it is the correct circle

                // if (std::abs((e->vertex(0)->coords() - m_O).norm() - m_R) < 1E-12)
                // {
                //     assert(std::abs((e->vertex(1)->coords() - m_O).norm() - m_R) < 1E-12);
                //     return false;
                // }
                    return false;
            }

            return true;
        }

    private:
        Eigen::Vector2d m_O;
        double m_R;
    };
}

//@}

#endif
