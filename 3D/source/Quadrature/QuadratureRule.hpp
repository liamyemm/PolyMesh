#include "GaussLegendre.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"
#include "Face.hpp"
#include "Cell.hpp"

#ifndef _QUADRATURERULE_HPP
#define _QUADRATURERULE_HPP

namespace PolyMesh3D
{
    namespace Quadrature
    {
        struct QuadratureNode
        {
            Eigen::Vector3d x;
            double w;
            QuadratureNode(Eigen::Vector3d x, double w) : x(x), w(w) {}
        };

        using QuadratureRule = std::vector<QuadratureNode>;

        const QuadratureRule generate_quadrature_rule(StraightMesh::Edge *edge, const GaussLegendre1D &quad);
        const QuadratureRule generate_quadrature_rule(StraightMesh::Face *face, const std::vector<QuadratureRule> &edge_quads, const GaussLegendre1D &quad);
        const QuadratureRule generate_quadrature_rule(StraightMesh::Cell *cell, const std::vector<QuadratureRule> &face_quads, const GaussLegendre1D &quad);
    }
}

#endif