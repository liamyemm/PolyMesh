#include "math.hpp"

namespace Math
{
    unsigned factorial(const unsigned n)
    {
        return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    }

    unsigned nChoosek(const unsigned n, unsigned k)
    {
        if (k > n)
        {
            return 0;
        }
        if (k * 2 > n)
        {
            k = n - k;
        }
        if (k == 0)
        {
            return 1;
        }

        size_t result = n;
        for (size_t i = 2; i <= k; ++i)
        {
            result *= (n - i + 1);
            result /= i;
        }
        return result;
    }

    // branch -1: -2\pi < \theta \le 0
    // branch 0:  -\pi < \theta \le \pi
    // branch 1:  0 \le \theta < 2\pi
    double atan2(const double y, const double x, int branch)
    {
        double val = (((x != 0) || (y != 0)) ? std::atan2(y, x) : 0.0);
        if (branch == 1)
        {
            if (y < 0)
            {
                val += 2.0 * PI;
            }
        }
        if (branch == -1)
        {
            if (y > 0)
            {
                val -= 2.0 * PI;
            }
            else if ((y == 0) && (x < 0))
            {
                val -= 2.0 * PI;
            }
        }
        return val;
    }

    // const double scalar_product(Eigen::MatrixXd A, Eigen::MatrixXd B)
    // {
    //     Eigen::Map<Eigen::VectorXd> vec_A(A.data(), A.cols() * A.rows());
    //     Eigen::Map<Eigen::VectorXd> vec_B(B.data(), B.cols() * B.rows());
    //     return vec_A.dot(vec_B);
    // }

    // const double scalar_product(Eigen::VectorXd A, Eigen::VectorXd B)
    // {
    //     return A.dot(B);
    // }

    const double scalar_product(double x, double y)
    {
        return x * y;
    }
    
    const double norm(double x)
    {
        return std::abs(x); 
    }
}
