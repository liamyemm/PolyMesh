#include "StokesTests.hpp"

namespace PolyMesh2D
{

    Eigen::Matrix2d tensor_product(const Eigen::Vector2d &v1, const Eigen::Vector2d &v2)
    {
        Eigen::Matrix2d tensorProduct;
        tensorProduct << v1(0) * v2(0), v1(0) * v2(1), v1(1) * v2(0), v1(1) * v2(1);
        return tensorProduct;
    }

    // ----------------------------------------------------------------------------
    //                          Implementation
    // ----------------------------------------------------------------------------

    StokesSingularity::StokesSingularity(double lambda, double alpha, double rot, Eigen::Vector2d x0) : m_lambda(lambda), m_alpha(alpha), m_rot(rot), m_x0(x0)
    {
    }

    Functional::Function<2, 1> StokesSingularity::invcurl_u() const
    {
        auto u_func = this->u();
        double lambda = m_lambda;

        std::function<double(Eigen::Vector2d)> value = [u_func, lambda](const Eigen::Vector2d &x) -> double
        {
            Eigen::Vector2d u = u_func.value(x);
            return 1.0 / (lambda + 1.0) * (x(1) * u(0) - x(0) * u(1));
        };
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> deriv = [u_func](const Eigen::Vector2d &x) -> Eigen::RowVector2d
        {
            Eigen::Vector2d u = u_func.value(x);
            return Eigen::RowVector2d(-u(1), u(0));
        };

        Functional::Function<2, 1> ret(Functional::Function<2, 1>(value, deriv));

        std::vector<Functional::Pole<Eigen::Vector2d>> poles;

        poles.push_back(Functional::Pole<Eigen::Vector2d>(m_x0, -1.0 - m_lambda));

        ret.set_poles(poles);

        return ret;
    }

    Functional::Function<2, 2> StokesSingularity::laplace_u() const
    {
        Eigen::Matrix2d rotationMatrix;
        rotationMatrix << std::cos(m_rot), -std::sin(m_rot), std::sin(m_rot), std::cos(m_rot);

        Eigen::Vector2d x0 = m_x0;
        double lambda = m_lambda;
        double alpha = m_alpha;

        std::function<Eigen::Vector2d(Eigen::Vector2d)> transform = [x0, rotationMatrix](const Eigen::Vector2d &x) -> Eigen::Vector2d
        {
            return rotationMatrix * (x - x0);
        };
        std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [rotationMatrix](const Eigen::Vector2d &x) -> Eigen::Matrix2d
        {
            return rotationMatrix;
        };

        std::function<Eigen::Vector2d(Eigen::Vector2d)> value = [lambda, alpha](const Eigen::Vector2d &x) -> Eigen::Vector2d
        {
            double r = x.norm();

            if (r < 1E-15)
            {
                if (lambda - 2 >= 0.0)
                {
                    return Eigen::RowVector2d::Zero();
                }
                else
                {
                    throw std::runtime_error("Attempting to evaluate function at singularity");
                }
            }

            double theta = Math::atan2(x(1), x(0), 1);

            double cos_l_alpha = std::cos(lambda * alpha);
            double sin_l_alpha = std::sin(lambda * alpha);
            double cos_alpha = std::cos(alpha);
            double sin_alpha = std::sin(alpha);

            double d1 = std::cos((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::sin((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha);

            double d2 = -std::sin((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::cos((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha);

            Eigen::Vector2d r_vec = x / r;
            Eigen::Vector2d theta_vec(-r_vec(1), r_vec(0));

            return (lambda - 1.0) * 2.0 * lambda * std::pow(r, lambda - 2) * (d1 * r_vec + d2 * theta_vec);
        };

        std::function<Eigen::Matrix2d(Eigen::Vector2d)> deriv = [lambda, alpha](const Eigen::Vector2d &x) -> Eigen::Matrix2d
        {
            double r = x.norm();

            if (r < 1E-15)
            {
                if (lambda - 3 >= 0.0)
                {
                    return Eigen::Matrix2d::Zero();
                }
                else
                {
                    throw std::runtime_error("Attempting to evaluate function at singularity");
                }
            }

            double theta = Math::atan2(x(1), x(0), 1);

            double cos_l_alpha = std::cos(lambda * alpha);
            double sin_l_alpha = std::sin(lambda * alpha);
            double cos_alpha = std::cos(alpha);
            double sin_alpha = std::sin(alpha);

            double d_lap_r_d_r = std::cos((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::sin((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha);

            double d_lap_theta_d_r = -std::sin((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::cos((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha);

            double d_lap_r_d_theta_minus_lap_theta = d_lap_theta_d_r;

            double d_lap_theta_d_theta_plus_lap_r = -d_lap_r_d_r;

            Eigen::Vector2d r_vec = x / r;
            Eigen::Vector2d theta_vec(-r_vec(1), r_vec(0));

            Eigen::Matrix2d r_r_mat = tensor_product(r_vec, r_vec);
            Eigen::Matrix2d r_theta_mat = tensor_product(r_vec, theta_vec);
            Eigen::Matrix2d theta_r_mat = tensor_product(theta_vec, r_vec);
            Eigen::Matrix2d theta_theta_mat = tensor_product(theta_vec, theta_vec);

            return (lambda - 2.0) * (lambda - 1.0) * 2.0 * lambda * std::pow(r, lambda - 3.0) * (d_lap_r_d_r * r_r_mat + d_lap_theta_d_r * theta_r_mat + d_lap_r_d_theta_minus_lap_theta * r_theta_mat + d_lap_theta_d_theta_plus_lap_r * theta_theta_mat);
        };

        Functional::Function<2, 2> ret(Functional::compose(Functional::Function<2, 2>(value, deriv), Functional::Function<2, 2>(transform, transform_deriv)));

        std::vector<Functional::Pole<Eigen::Vector2d>> poles;

        poles.push_back(Functional::Pole<Eigen::Vector2d>(x0, 2.0 - lambda));

        ret.set_poles(poles);

        return ret;
    }

    Functional::Function<2, 2> StokesSingularity::u() const
    {
        Eigen::Matrix2d rotationMatrix;
        rotationMatrix << std::cos(m_rot), -std::sin(m_rot), std::sin(m_rot), std::cos(m_rot);

        Eigen::Vector2d x0 = m_x0;
        double lambda = m_lambda;
        double alpha = m_alpha;

        std::function<Eigen::Vector2d(Eigen::Vector2d)> transform = [x0, rotationMatrix](const Eigen::Vector2d &x) -> Eigen::Vector2d
        {
            return rotationMatrix * (x - x0);
        };
        std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [rotationMatrix](const Eigen::Vector2d &x) -> Eigen::Matrix2d
        {
            return rotationMatrix;
        };

        std::function<Eigen::Vector2d(Eigen::Vector2d)> value = [lambda, alpha](const Eigen::Vector2d &x) -> Eigen::Vector2d
        {
            double r = x.norm();

            if (r < 1E-15)
            {
                if (lambda >= 0.0)
                {
                    return Eigen::Vector2d::Zero();
                }
                else
                {
                    throw std::runtime_error("Attempting to evaluate function at singularity");
                }
            }

            double theta = Math::atan2(x(1), x(0), 1);

            Eigen::Vector2d r_vec = x / r;
            Eigen::Vector2d theta_vec(-r_vec(1), r_vec(0));

            double cos_l_theta = std::cos(lambda * theta);
            double sin_l_theta = std::sin(lambda * theta);
            double cos_theta = std::cos(theta);
            double sin_theta = std::sin(theta);

            double cos_l_alpha = std::cos(lambda * alpha);
            double sin_l_alpha = std::sin(lambda * alpha);
            double cos_alpha = std::cos(alpha);
            double sin_alpha = std::sin(alpha);

            double tau_theta = (lambda * cos_l_theta * sin_theta - sin_l_theta * cos_theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - (sin_l_theta * sin_theta) / (sin_l_alpha * sin_alpha);

            double tau_r = (lambda * cos_l_theta * sin_theta + sin_l_theta * cos_theta) / ((lambda + 1) * sin_l_alpha * sin_alpha) + ((lambda - 1) * sin_l_theta * sin_theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha);

            return std::pow(r, lambda) * (tau_r * r_vec + tau_theta * theta_vec);
        };

        std::function<Eigen::Matrix2d(Eigen::Vector2d)> deriv = [lambda, alpha](const Eigen::Vector2d &x) -> Eigen::Matrix2d
        {
            double r = x.norm();

            if (r < 1E-15)
            {
                if (lambda - 1 >= 0.0)
                {
                    return Eigen::Matrix2d::Zero();
                }
                else
                {
                    throw std::runtime_error("Attempting to evaluate function at singularity");
                }
            }

            double theta = Math::atan2(x(1), x(0), 1);

            double cos_l_theta = std::cos(lambda * theta);
            double sin_l_theta = std::sin(lambda * theta);
            double cos_theta = std::cos(theta);
            double sin_theta = std::sin(theta);

            double cos_l_alpha = std::cos(lambda * alpha);
            double sin_l_alpha = std::sin(lambda * alpha);
            double cos_alpha = std::cos(alpha);
            double sin_alpha = std::sin(alpha);

            double d_tau_r_d_r = (lambda * cos_l_theta * sin_theta + sin_l_theta * cos_theta) / ((lambda + 1) * sin_l_alpha * sin_alpha) + ((lambda - 1) * sin_l_theta * sin_theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha);

            double d_tau_theta_d_r = (lambda * cos_l_theta * sin_theta - sin_l_theta * cos_theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - (sin_l_theta * sin_theta) / (sin_l_alpha * sin_alpha);

            double d_tau_r_d_theta_minus_tau_theta = d_tau_theta_d_r + 2.0 * ((cos_l_theta * cos_theta + sin_l_theta * sin_theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha) - (cos_l_theta * sin_theta - sin_l_theta * cos_theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha));

            double d_tau_theta_d_theta_plus_tau_r = -d_tau_r_d_r;

            Eigen::Vector2d r_vec = x / r;
            Eigen::Vector2d theta_vec(-r_vec(1), r_vec(0));

            Eigen::Matrix2d r_r_mat = tensor_product(r_vec, r_vec);
            Eigen::Matrix2d r_theta_mat = tensor_product(r_vec, theta_vec);
            Eigen::Matrix2d theta_r_mat = tensor_product(theta_vec, r_vec);
            Eigen::Matrix2d theta_theta_mat = tensor_product(theta_vec, theta_vec);

            return lambda * std::pow(r, lambda - 1) * (d_tau_r_d_r * r_r_mat + d_tau_r_d_theta_minus_tau_theta * r_theta_mat + d_tau_theta_d_r * theta_r_mat + d_tau_theta_d_theta_plus_tau_r * theta_theta_mat);
        };

        Functional::Function<2, 2> ret(Functional::compose(Functional::Function<2, 2>(value, deriv), Functional::Function<2, 2>(transform, transform_deriv)));

        std::vector<Functional::Pole<Eigen::Vector2d>> poles;

        poles.push_back(Functional::Pole<Eigen::Vector2d>(x0, -lambda));

        ret.set_poles(poles);

        return ret;
    }

    Functional::Function<2, 1> StokesSingularity::p(double avg_pressure) const
    {
        Eigen::Matrix2d rotationMatrix;
        rotationMatrix << std::cos(m_rot), -std::sin(m_rot), std::sin(m_rot), std::cos(m_rot);

        Eigen::Vector2d x0 = m_x0;
        double lambda = m_lambda;
        double alpha = m_alpha;

        std::function<Eigen::Vector2d(Eigen::Vector2d)> transform = [x0, rotationMatrix](const Eigen::Vector2d &x) -> Eigen::Vector2d
        {
            return rotationMatrix * (x - x0);
        };
        std::function<Eigen::Matrix2d(Eigen::Vector2d)> transform_deriv = [rotationMatrix](const Eigen::Vector2d &x) -> Eigen::Matrix2d
        {
            return rotationMatrix;
        };

        std::function<double(Eigen::Vector2d)> value = [lambda, alpha, avg_pressure](const Eigen::Vector2d &x) -> double
        {
            double r = x.norm();

            if (r < 1E-15)
            {
                if (lambda - 1 >= 0.0)
                {
                    return 0.0;
                }
                else
                {
                    throw std::runtime_error("Attempting to evaluate function at singularity");
                }
            }

            double theta = Math::atan2(x(1), x(0), 1);

            double cos_l_alpha = std::cos(lambda * alpha);
            double sin_l_alpha = std::sin(lambda * alpha);
            double cos_alpha = std::cos(alpha);
            double sin_alpha = std::sin(alpha);

            return 2.0 * lambda * std::pow(r, lambda - 1.0) * (std::cos((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::sin((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha)) - avg_pressure;
        };

        std::function<Eigen::RowVector2d(Eigen::Vector2d)> deriv = [lambda, alpha](const Eigen::Vector2d &x) -> Eigen::RowVector2d
        {
            double r = x.norm();

            if (r < 1E-15)
            {
                if (lambda - 2 >= 0.0)
                {
                    return Eigen::RowVector2d::Zero();
                }
                else
                {
                    throw std::runtime_error("Attempting to evaluate function at singularity");
                }
            }

            double theta = Math::atan2(x(1), x(0), 1);

            double cos_l_alpha = std::cos(lambda * alpha);
            double sin_l_alpha = std::sin(lambda * alpha);
            double cos_alpha = std::cos(alpha);
            double sin_alpha = std::sin(alpha);

            double d1 = std::cos((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::sin((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha);

            double d2 = -std::sin((lambda - 1.0) * theta) / (lambda * cos_l_alpha * sin_alpha - sin_l_alpha * cos_alpha) - std::cos((lambda - 1.0) * theta) / ((lambda + 1.0) * sin_l_alpha * sin_alpha);

            Eigen::RowVector2d r_vec = x.transpose() / r;
            Eigen::RowVector2d theta_vec(-r_vec(1), r_vec(0));

            return (lambda - 1.0) * 2.0 * lambda * std::pow(r, lambda - 2.0) * (d1 * r_vec + d2 * theta_vec);
        };

        Functional::Function<2, 1> ret(Functional::compose(Functional::Function<2, 1>(value, deriv), Functional::Function<2, 2>(transform, transform_deriv)));

        ret.add_pole(Functional::Pole<Eigen::Vector2d>(x0, 1.0 - lambda));

        return ret;
    }
}