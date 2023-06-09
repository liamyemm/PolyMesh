// for (int k = 0; k < n; k++)
// {
//     double old_z = 0.0;
//     double z = 2.0 * std::cos(Math::PI * (4.0 * k - 1.0 + 2.0 * alpha) / (4.0 * n + 2.0 * (alpha + beta + 1.0))) - 1.0; // initial guess for root k
//     double deriv_shifted_pn_z; // declared outside Newton loop, so we can use its final value in calculation of weights
//     for (its = 1; its <= MAXIT; its++)
//     {
//         // Newton method
//         if (std::abs(z - old_z) <= EPS)
//         {
//             break;
//         }

//         // evaluate pn and pn'
//         double shifted_pn_z = 0.0;
//         deriv_shifted_pn_z = 0.0;
//         for (int s = 0; s <= n; ++s)
//         {
//             double coeff = std::exp(std::lgamma(1 + alpha + n) + std::lgamma(1 + beta + n) - std::lgamma(1 + n - s) - std::lgamma(1 + beta + n - s) - std::lgamma(1 + s) - std::lgamma(1 + alpha + s));
//             shifted_pn_z += coeff * std::pow(z - 1, s) * std::pow(z, n - s);
//             if (s <= n - 1)
//             {
//                 double coeff_deriv = (1.0 + alpha + beta + n) * std::exp(std::lgamma(1 + alpha + n) + std::lgamma(1 + beta + n) - std::lgamma(n - s) - std::lgamma(1 + beta + n - s) - std::lgamma(1 + s) - std::lgamma(2 + alpha + s));
//                 deriv_shifted_pn_z += coeff_deriv * std::pow(z - 1, s) * std::pow(z, n - 1 - s);
//             }
//         }

//         old_z = z;
//         z = old_z - shifted_pn_z / deriv_shifted_pn_z; // Newton iteration
//     }

//     if (its > MAXIT)
//     {
//         throw("too many iterations in gauss_jacobi");
//     }
//     x[k] = z; // Store the root and the weight.
// }

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

QuadratureRule<double> gauss_jacobi(double alpha, double beta, int n)
{
    Eigen::MatrixXd comrade_mat = Eigen::MatrixXd::Zero(n, n); // maybe consider a sparse matrix as comrade_mat is tri-diagonal // <-- That is excessive. This method will probably be called with n around 5 or 10 and rarely more than 20.

    for (int k = 0; k < n - 1; ++k)
    {
        comrade_mat(k, k) = jacobi_comrade_diag(alpha, beta, k + 1);
        comrade_mat(k, k + 1) = comrade_mat(k + 1, k) = jacobi_comrade_subdiag(alpha, beta, k + 1);
    }
    comrade_mat(n - 1, n - 1) = jacobi_comrade_diag(alpha, beta, n);

    // Compute eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(comrade_mat);

    QuadratureRule<double> quad;
    quad.points = Eigen::Map<Eigen::VectorXd>(eig_solver.eigenvalues().data(), n);
    quad.weights = Eigen::Map<Eigen::ArrayXd>(eig_solver.eigenvectors().row(0).array().square().data(), n); // inefficient as only first entry of each eigen vector is used

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