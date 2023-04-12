#include "function.hpp"
#include "basis.hpp"
#include <array>
#include <Eigen/Dense>

#ifndef _COMMON_HPP
#define _COMMON_HPP

namespace PolyMesh2D
{
    namespace Functional
    {
        /*!
         * \addtogroup Functional
         * @{
         */

        /** 
         * Constructs a Curve object describing an arclength parameterisation of a line segment starting at v0 and terminating at v1.
        **/
        const Curve StraightLine(const ColVector &v0, const ColVector &v1); 

        /** 
         * Constructs a scalar-valued monomial function in 2 dimensions. 
         * @param powers is an array of two positive integer values storing the powers of each of the components of the monomial.
        **/
        const ScalarFunction2D Monomial(const std::array<size_t, 2> &powers);

        /** 
         * Constructs a scalar-valued monomial function in 1 dimension. 
         * @param power is a positive integer value which determines the power of the monomial.
        **/
        const ScalarFunction1D Monomial(const size_t power);

        /** 
         * Constructs a vector-valued monomial function in 2 dimensions. 
         * Only one component of the vector is non zero. @param index determines which vector component is non-zero
         * @param powers is an array of two positive integer values storing the powers of each of the components of the monomial.
        **/
        const VectorFunction2D VectorMonomial(const std::array<size_t, 2> &powers, const size_t index);     

        /** 
         * Generates a scalar-valued monomial basis in 2 dimensions.
         * @param transform describes a coordinate transformation under which the monomial is described. i.e. if the transformation is \f$\mathbf{y} = \mathbf{f}(\mathbf{x})\f$, then each basis function is a monomial in \f$\mathbf{y}\f$. transform is required to be affine (and non-degenerate) for the basis to span the same space as the polynomials in \f$\mathbf{x}\f$. No checks are performed to verify this.
         * @param degree determines the degree of the polynomial space.
        **/
        const ScalarBasis2D MonomialScalarBasis(const VectorFunction2D &transform, const size_t degree);

        /** 
         * Generates a vector-valued monomial basis in 2 dimensions.
         * @param transform describes a coordinate transformation under which the monomial is described. i.e. if the transformation is \f$\mathbf{y} = \mathbf{f}(\mathbf{x})\f$, then each basis function is a monomial in \f$\mathbf{y}\f$. transform is required to be affine (and non-degenerate) for the basis to span the same space as the polynomials in \f$\mathbf{x}\f$. No checks are performed to verify this.
         * @param degree determines the degree of the polynomial space.
        **/
        const VectorBasis2D MonomialVectorBasis(const VectorFunction2D &transform, const size_t degree);

        /** 
         * Generates a basis for the Raviart-Thomas-Nedelec space \f$\mathbb{RTN}^{k}(T) = \mathbb{P}^{k-1}(T)^d + \mathbf{x}\mathbb{P}^{k-1}(T)\f$ in 2 dimensions.
         * @param transform describes a coordinate transformation under which the basis is described. i.e. if the transformation is \f$\mathbf{y} = \mathbf{f}(\mathbf{x})\f$, then the basis for the the Raviart-Thomas-Nedelec is given by vector valued monomials in \f$\mathbf{y}\f$ plus \f$\mathbf{y}\f$ times scalar valued monomials in \f$\mathbf{y}\f$. transform is required to be affine (and non-degenerate) for the basis to span the same space as the Raviart-Thomas-Nedelec space in \f$\mathbf{x}\f$. No checks are performed to verify this.
         * @param degree determines the degree of the Raviart-Thomas-Nedelec space.
        **/
        const VectorBasis2D RaviartThomasNedelecBasis(const VectorFunction2D &transform, const size_t degree);

        /** 
         * Generates a set of functions which span \f$\mathbb{P}^0(E) + \mathbb{P}^k(\mathbb{R}^2)^2\cdot\mathbf{n}_E \f$ on a given Curve \f$E\f$.
         * @param param The parameterisation of the Curve \f$E\f$. 
         * @param degree The degree of the polynomial space.
        **/
        const ScalarBasis1D CurvedEdgeBasis(const Curve &param, const size_t degree);

        /** 
         * Generates a scalar-valued shifted monomial basis in the parameter of a Curve. i.e. If the Curve is defined by \f$E=\{\gamma(t)\}\f$ then the basis is in the coordinate \f$t\f$.
         * @param param The parameterisation of the Curve \f$E\f$. Required to be affine (i.e. a straight line segment) in order for the space to be equivalent to the restriction of the polynomial space \f$\mathbb{P}^k(\mathbb{R}^2)\f$ to the Curve \f$E\f$.
         * @param degree The degree of the polynomial space.
        **/
        const ScalarBasis1D MonomialScalarBasis(const Curve &param, const size_t degree);

        /** 
         * Constructs a Curve which is the restriction to a sub-domain of another Curve.
         * @param curve The original Curve.
         * @param t0 The new tmin for the constructed Curve.
         * @param t1 The new tmax for the constructed Curve.
        **/
        const Curve restriction(const Curve &curve, const double t0, const double t1);

        /** 
         * Uses a bisection method to find a point on a Curve which minimises the distance to a given point. This method assumes the signed curvature does not change sign so that
         * there exists at most one local maximum and one local minimum (given an infintite domain extension of the curve)
         * @param curve The Curve we are minimising the distance to.
         * @param point The point we are minimising the distance to.
         * @return The parameter value of the curve which minimises the distance to point.
        **/
        const double min_distance(const Curve &curve, const Eigen::Vector2d &point);

        //@}
    }
}
#endif
