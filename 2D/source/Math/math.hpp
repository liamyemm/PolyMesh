#ifndef _MATH_HPP
#define _MATH_HPP

#include <cmath>
#include <Eigen/Dense>

namespace Math
{
    /// Free math functions and global variables ///

    const double PI = acos(-1);

    unsigned factorial(const unsigned n);

    unsigned nChoosek(const unsigned n, unsigned k);

    template <typename T>
    int sgn(T val);

    // branch -1: -2\pi < \theta \le 0
    // branch 0:  -\pi < \theta \le \pi
    // branch 1:  0 \le \theta < 2\pi
    double atan2(const double y, const double x, int branch = 0);

    const double scalar_product(double x, double y);

    template <int n_rows, int n_cols>
    const double scalar_product(const Eigen::Matrix<double, n_rows, n_cols> &x, const Eigen::Matrix<double, n_rows, n_cols> &y);

    const double norm(double x);

    template <int n_rows, int n_cols>
    const double norm(const Eigen::Matrix<double, n_rows, n_cols> &x);
}

// Implementation of templated function

template <typename T>
int Math::sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <int n_rows, int n_cols>
const double Math::scalar_product(const Eigen::Matrix<double, n_rows, n_cols> &x, const Eigen::Matrix<double, n_rows, n_cols> &y)
{
    double val = 0.0;
    for (int i = 0; i < n_rows; ++i)
    {
        for (int j = 0; j < n_cols; ++j)
        {
            val += x(i, j) * y(i, j);
        }
    }
    return val;
}

template <int n_rows, int n_cols>
const double Math::norm(const Eigen::Matrix<double, n_rows, n_cols> &x)
{
    return std::sqrt(Math::scalar_product(x, x));
}

#endif
