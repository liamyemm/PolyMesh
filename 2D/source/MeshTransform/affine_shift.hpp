#include <Eigen/Core>

#ifndef _AFFINE_SHIFT_HPP
#define _AFFINE_SHIFT_HPP

namespace PolyMesh2D
{
    // namespace MeshTransform
    // {
        void affine_shift(const std::string &input_file_path, const std::string &output_file_path, const Eigen::Matrix2d &A, const Eigen::Vector2d &b);
    // }
}

#endif