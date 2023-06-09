#include "function.hpp"
#include "basis.hpp"
#include "QuadratureRule.hpp"
#include "generate_quadrature_rule.hpp"
#include "QuadHandler.hpp"

#include <Mesh>
#include <Eigen/Dense>

#ifndef _INTEGRATION_HPP
#define _INTEGRATION_HPP

namespace PolyMesh2D
{
    namespace Quadrature
    {
        template<unsigned output_dim>
        Eigen::MatrixXd integrate(const Basis<2, output_dim> &basis1, const Basis<2, output_dim> &basis2, const Cell * cell, const QuadHandler &quad)
        {
            Eigen::MatrixXd mat(basis1.dimension(), basis2.dimension());

            for(size_t i = 0; i < basis1.dimension(); ++i)
            {
                for(size_t j = 0; j < basis2.dimension(); ++j)
                {
                    mat(i, j) = quad.integrate( PolyMesh2D::functional::scalar_product(basis1.function(i), basis2.function(j)) );
                }
            }

            return mat;
        }
    }
}

#endif