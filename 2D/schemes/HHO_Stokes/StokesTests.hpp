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
}

//@}

#endif
