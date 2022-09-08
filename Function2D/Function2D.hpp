// standard libraries
#include <vector>    // std::vector
#include <algorithm> // std::find
#include <fstream> // std::ofstream
 
// external linking to Eigen is required
#include <Eigen/Dense> // Eigen::Matrix, Eigen::MatrixXd, Eigen::FullPivLU

// path to Math specified so that external linking is NOT required
#include <../Math/math.hpp> // Math::factorial, Math::sgn

#include <../Mesh2D/Polytope2D.hpp>

#ifndef _FUNCTION2D_HPP
#define _FUNCTION2D_HPP

namespace Function2D
{
    enum xy
    {
        x = 0,
        y = 1
    };

    class ScalarFunction1D // maps scalar to scalar
    {
        public:
            virtual double value(double t) = 0;
            virtual double derivative(double t) = 0;
    };

    class ScalarFunction2D // maps vector to scalar
    {
        public:
            virtual double value(Mesh2D::VectorRd x) = 0;
            virtual Mesh2D::VectorRd gradient(Mesh2D::VectorRd x) = 0;
    };

    class VectorFunction // maps vector to vector
    {
        public:
            virtual Mesh2D::VectorRd value(Mesh2D::VectorRd x) = 0;
            virtual Mesh2D::MatrixRd gradient(Mesh2D::VectorRd x) = 0;
    };

    class ParametricFunction // maps scalar to vector
    {
        public:
            virtual Mesh2D::VectorRd value(double t) = 0;
            virtual Mesh2D::VectorRd derivative(double t) = 0;

            const double tmin;
            const double tmax;

        protected:
            bool input_check(double t)
            {
                return (tmin <= t && t <= tmax);
            }
    };

    // class InverseDivergence : VectorFunction // given a scalar function u, returns the function a_r e_r with e_r the unit radial vector and a_r = r^{-1} \int_{0}^1 s \hat{u}(s, \theta) ds and \hat{u}(r, \theta) = u(r \cos (\theta), r \sin (\theta))
    // {
    //     public:
    //         InverseDivergence(ScalarFunction2D *scalar_func) : _scalar_func(scalar_func) {}

    //         Mesh2D::VectorRd value(Mesh2D::VectorRd x);
    //         Mesh2D::MatrixRd gradient(Mesh2D::VectorRd x);
    //     private:
    //         ScalarFunction2D *_scalar_func
    // }

    class StraightEdgeParametrisation : ParametricFunction
    {
        public:
            StraightEdgeParametrisation(Mesh2D::Edge * edge) : v0(edge->get_simplices()[0][0]) , v1(edge->get_simplices()[0][1]) , tmin(0.0) , tmax(1.0) {}

            Mesh2D::VectorRd value(double t) 
            {
                assert(input_check(t));
                return v0 + t * (v1 - v0);
            }

            Mesh2D::VectorRd derivative(double t) 
            {
                assert(input_check(t));
                return v1 - v0;
            }

        private:
            Mesh2D::VectorRd v0;
            Mesh2D::VectorRd v1;
    };

    // class Monomial : ScalarFunction2D
    // {
    //     public:
    //         Monomial(int degree, xy coord) : _degree(degree) , _coord(coord) 
    //         {
    //             assert(_degree >= 0);
    //         }

    //         double value(Mesh2D::VectorRd x)
    //         {
    //             return std::pow(x(_coord), _degree);
    //         }

    //         Mesh2D::VectorRd gradient(Mesh2D::VectorRd x)
    //         {
    //             VectorRd grad(0.0, 0.0);
    //             grad(coord) = (_degree == 0 ? 0.0 : double(_degree) * std::pow(x(_coord), _degree - 1));
    //             return grad;
    //         }

    //     private: 
    //         int _degree;
    //         xy _coord;
    // };

    inline double diameter(ParametricFunction *func, int partitions)
    {
        double deltat = (func->tmax - func->tmin) / double(partitions);
        double diam = 0.0;

        for(int i = 0; i < partitions; ++i)
        {
            for(int j = i + 1; j < partitions; ++j)
            {
                double ti = func->tmin + double(i) * deltat;
                double tj = func->tmin + double(j) * deltat;
                diam = std::max(diam, (func->value(ti) - func->value(tj)).norm());
            }
        }

        return diam;
    }
}


#endif
