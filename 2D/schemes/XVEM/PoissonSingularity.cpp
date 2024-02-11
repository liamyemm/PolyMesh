#include "PoissonSingularity.h"

namespace PolyMesh2D
{
    Eigen::Matrix2d rotationMatrix(double phi)
    {
        Eigen::Matrix2d mat;

        mat << std::cos(phi), -std::sin(phi),
               std::sin(phi),  std::cos(phi);

        return mat;
    }

    const Functional::Function<2, 1> PoissonSingularity(double alpha, double rot, const Eigen::Vector2d &shift)
    {
        // std::cout << "\n" << rot << "\n";
        // std::cout << "\n\n" << shift << "\n\n";
        Eigen::Matrix2d rotMat = rotationMatrix(rot);

        std::function<double(Eigen::Vector2d)> value = [alpha, rotMat, shift](const Eigen::Vector2d &x) -> double
        {
            Eigen::Vector2d y = rotMat * x - shift;
            double r = y.norm();
            double theta = Math::atan2(y(1), y(0), 1);

            return std::pow(r, alpha) * std::sin(alpha * theta);
        };
        
        std::function<Eigen::RowVector2d(Eigen::Vector2d)> deriv = [alpha, rotMat, shift](const Eigen::Vector2d &x) -> Eigen::RowVector2d
        {
            Eigen::Vector2d y = rotMat * x - shift;
            double r = y.norm();
            double theta = Math::atan2(y(1), y(0), 1);

            return alpha * std::pow(r, alpha - 1.0) * Eigen::RowVector2d(std::sin((alpha - 1.0) * theta), std::cos((alpha - 1.0) * theta)) * rotMat; 
        };

        Functional::Function<2, 1> ret(Functional::Function<2, 1>(value, deriv));

        ret.add_pole(Functional::Pole<Eigen::Vector2d>(shift, -alpha));

        return ret;
    }
}