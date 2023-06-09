#include "function.hpp"

namespace PolyMesh2D
{
    namespace Functional
    {
        Curve::Curve(const double t0, const double t1, const std::function<Eigen::Vector2d(double)> &val, const std::function<Eigen::Vector2d(double)> &deriv) : VectorFunction1D(val, deriv), tmin(t0), tmax(t1)
        {
            // __value = val;
            // __derivative = deriv;
        }

        ColVector Curve::value(const double &t) const
        {
            assert(input_check(t));
            return VectorFunction1D::value(t);
        }
        ColVector Curve::derivative(const double &t) const
        {
            assert(input_check(t));
            return VectorFunction1D::derivative(t);
        }

        ColVector Curve::tangent(const double t) const
        {
            return this->derivative(t).normalized();
        }

        ColVector Curve::normal(const double t) const
        {
            ColVector tan(this->tangent(t));
            return ColVector(-tan(1), tan(0));
        }

        double Curve::derivative_norm(const double t) const
        {
            return this->derivative(t).norm();
        }

        // double Curve::curvature(const double t) const
        // {
        //     assert(input_check(t));
        //     return __curvature(t);
        // }

        bool Curve::input_check(const double t) const
        {
            return (tmin - 1E-15 < t && t < tmax + 1E-15);
        }
    }
}