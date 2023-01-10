#ifndef _MATH_HPP
#define _MATH_HPP

#include <cmath>

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
}

// Implementation of templated function

template <typename T>
int Math::sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

#endif
