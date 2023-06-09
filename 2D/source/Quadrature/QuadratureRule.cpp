#include "QuadratureRule.hpp"

namespace PolyMesh2D
{
    namespace Quadrature
    {
        // https://www.ams.org/journals/mcom/1969-23-106/S0025-5718-69-99647-1/S0025-5718-69-99647-1.pdf
        
        // method to return the kth subdiagonal entry of the comrade matrix associated with Jacobi polynomials
        double jacobi_comrade_subdiag(double alpha, double beta, int k)
        {
            assert(k >= 1);
            assert(alpha > -1);
            assert(beta > -1);

            double a = alpha + k;
            double b = beta + k;
            double c = a + b;
            if (k == 1)
            {
                return std::sqrt(a * b / (c + 1.0)) / c; // avoids division by zero at alpha + beta = -1 (c = 1)
            }
            return std::sqrt(k * a * b * (c - k) / (c * c - 1.0)) / c;
        }

        // method to return the kth diagonal entry of the comrade matrix associated with Jacobi polynomials
        double jacobi_comrade_diag(double alpha, double beta, int k)
        {
            assert(k >= 1);
            assert(alpha > -1);
            assert(beta > -1);

            if (k == 1)
            {
                return 0.5 * (1.0 + (beta - alpha) / (2.0 + alpha + beta)); // avoids division by zero at alpha + beta = 0
            }
            return 0.5 * (1.0 + (beta - alpha) * (beta + alpha) / ((2 * k + alpha + beta) * (2 * k + alpha + beta - 2.0)));
        }

        const QuadratureRule<double> gauss_jacobi(const double alpha, const double beta, const unsigned doe)
        {
            unsigned n = doe % 2 == 0 ? (doe + 2) / 2 : (doe + 1) / 2; // std::ceil((doe + 1) / 2.0)
            Eigen::MatrixXd comrade_mat = Eigen::MatrixXd::Zero(n, n); // maybe consider a sparse matrix as comrade_mat is tri-diagonal // <-- That is excessive. This method will probably be called with n around 5 or 10 and rarely more than 20.

            for (unsigned k = 0; k < n - 1; ++k)
            {
                comrade_mat(k, k) = jacobi_comrade_diag(alpha, beta, k + 1);
                comrade_mat(k, k + 1) = comrade_mat(k + 1, k) = jacobi_comrade_subdiag(alpha, beta, k + 1);
            }
            comrade_mat(n - 1, n - 1) = jacobi_comrade_diag(alpha, beta, n);

            // Compute eigenvalues
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(comrade_mat);

            Eigen::VectorXd points = eig_solver.eigenvalues();
            Eigen::ArrayXd weights = std::exp(std::lgamma(alpha + 1) + std::lgamma(beta + 1) - std::lgamma(alpha + beta + 2)) * (eig_solver.eigenvectors().row(0).array().square());

            QuadratureRule<double> quad;
            quad.points.reserve(n);
            quad.weights.reserve(n);

            for (unsigned i = 0; i < n; ++i)
            {
                quad.points.push_back(points(i));
                quad.weights.push_back(weights(i));
            }

            return quad;
        }

        void normalise_weights(QuadratureRule<double> &quad, double alpha, double beta)
        {
            assert(alpha > -1);
            assert(beta > -1);

            for (size_t i = 0; i < quad.size(); ++i)
            {
                double xi = quad.point(i);
                double denominator = std::pow(xi, beta) * std::pow(1.0 - xi, alpha);
                quad.weights[i] /= denominator;
            }
        }
    }
}
