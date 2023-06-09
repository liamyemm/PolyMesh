#include <Eigen/Dense>

#include "math.hpp"

#ifndef _QUADRATURERULE_HPP
#define _QUADRATURERULE_HPP

namespace PolyMesh2D
{
    namespace Quadrature
    {
        /**
         * @brief A structure representing a quadrature node, consisting of a point and a weight.
         * @tparam PointType The type of the node point, which can be doubles or Eigen::Vector2d.
         */
        template <typename PointType>
        struct QuadratureNode
        {
            PointType x; /**< The node point. */
            double w;    /**< The node weight. */

            /**
             * @brief Constructs a new QuadratureNode object.
             * @param x The node point.
             * @param w The node weight.
             */
            QuadratureNode(PointType x, double w) : x(x), w(w) {}
        };

        /**
         * @brief A structure representing a quadrature rule, consisting of a set of nodes (points and weights).
         * @tparam PointType The type of the node points, which can be doubles or Eigen::Vector2d.
         */
        template <typename PointType>
        struct QuadratureRule
        {
            std::vector<PointType> points; /**< The set of node points. */
            std::vector<double> weights;   /**< The set of node weights. */

            /**
             * @brief Returns the i-th node point in the quadrature rule.
             * @param i The index of the node point to retrieve.
             * @throws std::out_of_range if i is out of range.
             * @return The i-th node point.
             */
            PointType point(int i) const
            {
                if (i < 0 || i >= static_cast<int>(points.size()))
                {
                    throw std::out_of_range("Index out of range");
                }
                return points[i];
            }

            /**
             * @brief Returns the weight of the i-th node in the quadrature rule.
             * @param i The index of the node weight to retrieve.
             * @throws std::out_of_range if i is out of range.
             * @return The weight of the i-th node.
             */
            double weight(int i) const
            {
                if (i < 0 || i >= static_cast<int>(weights.size()))
                {
                    throw std::out_of_range("Index out of range");
                }
                return weights[i];
            }

            /**
             * @brief Returns the number of nodes in the quadrature rule.
             * @return The number of nodes in the quadrature rule.
             */
            size_t size() const
            {
                return points.size();
            }

            /**
             * @brief Returns the i-th node in the quadrature rule as a QuadratureNode.
             * @param i The index of the node to retrieve.
             * @throws std::out_of_range if i is out of range.
             * @return The i-th node as a QuadratureNode.
             * @note This function is included for compatibility with old code using QuadratureRule = std::vector<QuadratureNode<PointType>>.
             */
            QuadratureNode<PointType> operator[](int i) const
            {
                if (i < 0 || i >= static_cast<int>(points.size()))
                {
                    throw std::out_of_range("Index out of range");
                }
                return QuadratureNode<PointType>(points[i], weights[i]);
            }
        };

        /**
         * @brief Computes a quadrature rule on (0, 1) using Gauss-Jacobi quadrature.
         *
         * @param alpha The paramater a of the Jacobi polynomial defined by the orthogonality relation \int_{0}^{1} x^b (1 - x)^a P_i^{a,b}(x) P_j^{a,b}(x) = 0.
         * @param beta The paramater b of the Jacobi polynomial defined by the orthogonality relation \int_{0}^{1} x^b (1 - x)^a P_i^{a,b}(x) P_j^{a,b}(x) = 0.
         * @param doe The degree of exactness of the quadrature rule.
         *
         * @return A quadrature rule consisting of n nodes and weights.
         */
        const QuadratureRule<double> gauss_jacobi(const double alpha, const double beta, const unsigned doe);

        /**
         * @brief Re-scales the weights of a Gauss-Jacobi rule with parameters alpha, beta.
         * 
         * The Gauss-Jacobi rule approximates integrals of the form \int_{0}^{1} x^b (1 - x)^a g(x) by the sum \sum_{i=1}^n w_i g(x_i). However, sometimes we want to write f(x) = x^b (1 - x)^a g(x) and evaluate f at the absiccae rather than g. For this, we require to re-scale each weight w_i by x_i^{-b} (1 - x_i)^{-a}.
         *
         * @param quad The Gauss-Jacobi rule to be re-scaled.
         * @param alpha The paramater a in the Jacobi polynomial P_n^{a,b}(x).
         * @param beta The paramater b in the Jacobi polynomial P_n^{a,b}(x).
         */
        void normalise_weights(QuadratureRule<double> &quad, double alpha, double beta);
    }
}

#endif