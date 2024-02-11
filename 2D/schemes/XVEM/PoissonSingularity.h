#ifndef _POISSONSINGULARITY_H
#define _POISSONSINGULARITY_H

#include "function.hpp"
#include "Mesh.hpp"
#include "math.hpp"

namespace PolyMesh2D
{
    const Functional::Function<2, 1> PoissonSingularity(double alpha, double rot = 0.0, const Eigen::Vector2d &shift = Eigen::Vector2d::Zero());

} // namespace PolyMesh2D

#endif // _POISSONSINGULARITY_H
