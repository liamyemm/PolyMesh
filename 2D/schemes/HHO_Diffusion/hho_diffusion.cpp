#include "hho_diffusion.hpp"

// internal libraries
#include "parallel_for.hpp"
#include "vtu_writer.hpp"
#include "QuadratureRule.hpp"
#include "basis.hpp"
#include "function.hpp"

#include "../HHO_Poisson/TestCase.hpp"

// mesh intersection
#include "../CutMesh/CutMesh.hpp"
#include "../IntersectMesh/IntersectMesh.hpp"

// boost libraries
#include <boost/program_options.hpp> // boost::program_options
#include <boost/timer/timer.hpp>     // boost::cpu_timer

// std libraries
#include <iostream> // std::cout
#include <iomanip>  // std::setw, std::left

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#endif

using namespace PolyMesh2D::HHODIFFUSION;

int main(const int argc, const char **argv)
{
    ModelParameters params(argc, argv);

    // Build mesh and reorder edges
    PolyMesh2D::StraightMesh::MeshBuilder builder = PolyMesh2D::StraightMesh::MeshBuilder(params.mesh_name);
    std::unique_ptr<PolyMesh2D::StraightMesh::Mesh> straight_mesh;

    try
    {
        straight_mesh = builder.build_the_mesh();
    }
    catch (std::string msg)
    {
        std::cerr << msg;
    }

    std::function<Eigen::Vector2d(double)> circle_val_outer = [](double t) -> Eigen::Vector2d
    { return Eigen::Vector2d(std::cos(t), std::sin(t)); };
    std::function<Eigen::Vector2d(double)> circle_deriv_outer = [](double t) -> Eigen::Vector2d
    { return Eigen::Vector2d(-std::sin(t), std::cos(t)); };

    PolyMesh2D::Functional::Curve bdry_param_outer(0.0, 2.0 * Math::PI, circle_val_outer, circle_deriv_outer);

    std::function<double(PolyMesh2D::Functional::ColVector)> LS_outer = [](const PolyMesh2D::Functional::ColVector &x) -> double
    {
        return 1.0 - (std::pow(x(0), 2) + std::pow(x(1), 2));
    };

    ScalarFunction2D level_set_outer(LS_outer);

    PolyMesh2D::MeshCutter mesh_cutter_outer(straight_mesh.get(), level_set_outer, bdry_param_outer);
    std::unique_ptr<Mesh> curved_mesh_outer = mesh_cutter_outer.cut_mesh();

    straight_mesh.reset(); // delete straight mesh

    double R = 0.8;

    std::function<Eigen::Vector2d(double)> circle_val = [&R](double t) -> Eigen::Vector2d
    { return R * Eigen::Vector2d(std::cos(t), std::sin(t)); };
    std::function<Eigen::Vector2d(double)> circle_deriv = [&R](double t) -> Eigen::Vector2d
    { return R * Eigen::Vector2d(-std::sin(t), std::cos(t)); };

    PolyMesh2D::Functional::Curve bdry_param(0.0, 2.0 * Math::PI, circle_val, circle_deriv);

    std::function<double(PolyMesh2D::Functional::ColVector)> LS = [&R](const PolyMesh2D::Functional::ColVector &x) -> double
    {
        return R * R - (std::pow(x(0), 2) + std::pow(x(1), 2));
    };

    ScalarFunction2D level_set(LS);

    // PolyMesh2D::MeshCutter2 mesh_cutter(straight_mesh.get(), level_set, bdry_param);
    PolyMesh2D::MeshCutter2 mesh_cutter(curved_mesh_outer.get(), level_set, bdry_param);
    std::unique_ptr<Mesh> curved_mesh = mesh_cutter.cut_mesh();

    // straight_mesh.reset(); // delete straight mesh

    // Reordering
    PolyMesh2D::HHOPOISSON::reorder_edges(curved_mesh.get());

    std::ofstream mesh_out("mesh_plot.dat");
    curved_mesh->plot_mesh(&mesh_out, 20);
    mesh_out.close();
    // exit(0);

    assert(curved_mesh->test());

    // Print the data
    std::cout << "[Scheme] Data:\n";
    std::cout << "     No. cells = " << curved_mesh->n_cells() << ", No. edges = " << curved_mesh->n_edges() << ", No. vertices = " << curved_mesh->n_vertices() << "\n";
    std::cout << "     Mesh = " << params.mesh_name << "\n";
    std::cout << "     Degrees: edge = " << params.edge_degree << "; cell = " << params.cell_degree << "\n";
    std::cout << "     Using threads = " << (params.use_threads ? "true" : "false") << std::endl;

    HybridCore hho(curved_mesh.get(), params.cell_degree, params.edge_degree, params.use_threads, params.orthonormalise);

    // std::function<double(Eigen::Vector2d, double)> src_helper = [](const Eigen::Vector2d &x, const double anisotropy) -> double
    // {
    //     double r_hat_sq = std::pow(x(0) + 0.25, 2) + std::pow(x(1) + 0.25, 2);
    //     double r_hat = std::sqrt(r_hat_sq);

    //     double sin_pi_x = std::sin(Math::PI * x(0));
    //     double sin_pi_y = std::sin(Math::PI * x(1));

    //     // std::cout << anisotropy << "\n";

    //     return (4.0 * Math::PI * (0.25 + x(1)) * (r_hat_sq) * (1.0 - r_hat) * std::cos(Math::PI * x(1)) * sin_pi_x - 2.0 * std::pow(0.25 + x(1), 2) * r_hat * sin_pi_x * sin_pi_y + 2.0 * (r_hat_sq) * (1.0 - r_hat) * sin_pi_x * sin_pi_y + std::pow(Math::PI, 2) * std::pow(r_hat, 3) * std::pow(1.0 - r_hat, 2) * sin_pi_x * sin_pi_y - 2.0 * std::pow(0.25 + x(1), 2) * (1.0 - r_hat) * sin_pi_x * sin_pi_y + anisotropy * (4 * Math::PI * (0.25 + x(0)) * (r_hat_sq) * (1.0 - r_hat) * std::cos(Math::PI * x(0)) - 2.0 * std::pow(0.25 + x(0), 2) * r_hat * sin_pi_x + 2.0 * (r_hat_sq) * (1.0 - r_hat) * sin_pi_x + std::pow(Math::PI, 2) * std::pow(r_hat, 3) * std::pow(1.0 - r_hat, 2) * sin_pi_x - 2.0 * std::pow(0.25 + x(0), 2) * (1.0 - r_hat) * sin_pi_x) * sin_pi_y) / std::pow(r_hat, 3);
    // };

    // std::function<double(Eigen::Vector2d, double)> src_helper = [](const Eigen::Vector2d &x, const double anisotropy) -> double
    // {
    //     return -(-24.0 * Math::PI * (1.0 + 4.0 * x(1)) * std::pow(-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2), 2) * std::cos(Math::PI * x(1)) * std::sin(Math::PI * x(0)) - 96.0 * std::pow(1.0 + 4.0 * x(1), 2) * (-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2)) * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)) - 48.0 * std::pow(-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2), 2) * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)) + std::pow(Math::PI, 2) * std::pow(-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2), 3) * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)) + 512.0 * anisotropy * (1.0 - std::pow(0.25 + x(0), 2) - std::pow(0.25 + x(1), 2)) * ((3.0 * Math::PI * (1.0 + 4.0 * x(0)) * (-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2)) * std::cos(Math::PI * x(0))) / 8.0 + 24.0 * std::pow(0.25 + x(0), 2) * std::sin(Math::PI * x(0)) + (3.0 * (-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2)) * std::sin(Math::PI * x(0))) / 4.0 - (std::pow(Math::PI, 2) * std::pow(-7.0 + 4.0 * x(0) + 8.0 * std::pow(x(0), 2) + 4.0 * x(1) + 8.0 * std::pow(x(1), 2), 2) * std::sin(Math::PI * x(0))) / 64.0) * std::sin(Math::PI * x(1))) / 512.0;
    // };

    // double anisotropy = 1E6;

    // double alpha1 = 1.0;
    // double alpha2 = 1E12;

    // double alpha = 1E3;
    double beta1 = 1E-6;
    double beta2 = 1.0;

    // std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [&src_helper, &anisotropy, &level_set](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    // {
    //     // double an = 1.0;

    //     double an;
    //     // if ((cell->center_mass() + Eigen::Vector2d(0.25, 0.25)).norm() < 1.0)
    //     // if ( std::abs(level_set.value(x)) < 1E-12)
    //     // {
    //     //     if (level_set.value(cell->center_mass()) < 0.0)
    //     //     {
    //     //         an = anisotropy;
    //     //     }
    //     //     else
    //     //     {
    //     //         an = 1.0;
    //     //     }
    //     //     // std::cout << "hello\n";
    //     // }
    //     // else
    //     if  (level_set.value(x) > 0.0)
    //     {
    //         an = anisotropy;
    //     }
    //     else
    //     {
    //         an = 1.0;
    //     }
    //     // else
    //     // {
    //     //     std::cout << "hello1\n";
    //     // }

    //     return src_helper(x, an);
    // };
    double a = 2.5;
    // std::function<double(Eigen::Vector2d)> u = [&beta1, &beta2, &R](const Eigen::Vector2d &x) -> double
    std::function<double(Eigen::Vector2d)> u = [&a, &R](const Eigen::Vector2d &x) -> double
    {
        double r_sq = x.squaredNorm();
        double R_sq = R * R;
        if (std::abs(r_sq - R_sq) < 1E-14)
        {
            return 0.0;
        }
        return (1.0 - r_sq) * std::pow(std::abs(R_sq - r_sq), a);
        // if (std::abs(x.norm() - R) < 1E-14)
        // {
        //     return 0.0;
        // }
        // double r_sq = x.squaredNorm();
        // double R_sq = R * R;

        // if (x.norm() < R)
        // {
        //     // r < R
        //     // return std::pow(R * R - x.squaredNorm(), 2) * (R * R - 4.0 * x.squaredNorm() - 6.0 * (1.0 - beta2) * x(0) * x(1));
        //     // return x(0) * x(1) * std::pow(r_sq - R_sq, 2) * (15 * x(0) * x(1) * (21.0 * std::pow(r_sq, 2) - 22.0 * r_sq * R_sq + 5.0 * std::pow(R_sq, 2)) + (1.0 - beta2) * (78.0 * std::pow(r_sq, 3) - 119.0 * std::pow(r_sq, 2) * R_sq + 44.0 * r_sq * std::pow(R_sq, 2) - 3.0 * std::pow(R_sq, 3) + 40.0 * std::pow(x(0), 2) * std::pow(x(1), 2) * (7.0 * r_sq - 4.0 * R_sq)));
        //     return std::pow(R_sq - r_sq, 2) * (1.0 - r_sq) * (1.0 + (r_sq / R_sq) * (6.0 * r_sq - (3.0 * R_sq + 4.0)) + 2.0 * (1.0 - beta2) * x(0) * x(1) * (5.0 * r_sq - (3.0 + 2.0 * R_sq)) / R_sq);
        // }
        // else
        // {
        //     // r > R
        //     // return std::pow(R * R - x.squaredNorm(), 2) * (R * R - 4.0 * x.squaredNorm() - 6.0 * (1.0 - beta1) * x(0) * x(1));
        //     return std::pow(R_sq - r_sq, 2) * (1.0 - r_sq) * (1.0 + (r_sq / R_sq) * (6.0 * r_sq - (3.0 * R_sq + 4.0)) + 2.0 * (1.0 - beta1) * x(0) * x(1) * (5.0 * r_sq - (3.0 + 2.0 * R_sq)) / R_sq);
        // }
    };
    std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Du = [&beta1, &beta2, &R](const Eigen::Vector2d &x) -> PolyMesh2D::Functional::RowVector
    {
        if (std::abs(x.norm() - R) < 1E-14)
        {
            return PolyMesh2D::Functional::RowVector(0.0, 0.0);
        }
        double r_sq = x.squaredNorm();
        double R_sq = R * R;
        if (x.norm() < R)
        {
            // r < R
            // return PolyMesh2D::Functional::RowVector(-6 * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * (-4 * std::pow(x(0), 3) + 5 * (-1 + beta2) * std::pow(x(0), 2) * x(1) - 4 * x(0) * std::pow(x(1), 2) + (-1 + beta2) * std::pow(x(1), 3) + std::pow(R, 2) * (2 * x(0) + x(1) - beta2 *, x(1))), 6 * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * (-((-1 + beta2) * std::pow(x(0), 3)) + std::pow(R, 2) * ((-1 + beta2) * x(0) - 2 * x(1)) + 4 * std::pow(x(0), 2) * x(1) - 5 * (-1 + beta2) * x(0) * std::pow(x(1), 2) + 4 * std::pow(x(1), 3)));

            // return (r_sq - R_sq) * PolyMesh2D::Functional::RowVector(x(1) * (3.0 * (1.0 - beta2) * std::pow(R, 8) + 50 * std::pow(R, 6) * x(0) * (2 * (-1.0 + beta2) * x(0) - 3 * x(1)) + 1200 * (-1.0 + beta2) * std::pow(R, 2) * std::pow(x(0), 4) * std::pow(x(1), 2) + 480 * std::pow(R, 4) * std::pow(x(0), 2) * x(1) * (2 * x(0) + x(1) - beta2 * x(1)) + (47 * (-1.0 + beta2) * std::pow(R, 6) - 1680 * (-1.0 + beta2) * std::pow(x(0), 4) * std::pow(x(1), 2) + 10 * std::pow(R, 4) * x(0) * (-74 * (-1.0 + beta2) * x(0) + 81 * x(1)) - 120 * std::pow(R, 2) * std::pow(x(0), 2) * x(1) * (27 * x(0) - 11 * (-1.0 + beta2) * x(1))) * r_sq + (-163 * (-1.0 + beta2) * std::pow(R, 4) + 10 * std::pow(R, 2) * x(0) * (142 * (-1.0 + beta2) * x(0) - 129 * x(1)) + 840 * std::pow(x(0), 2) * x(1) * (3 * x(0) + x(1) - beta2 * x(1))) * std::pow(r_sq, 2) + (197 * (-1.0 + beta2) * std::pow(R, 2) + 30 * x(0) * (-26 * (-1.0 + beta2) * x(0) + 21 * x(1))) * std::pow(r_sq, 3) - 78 * (-1.0 + beta2) * std::pow(r_sq, 4)), x(0) * (-3 * (-1.0 + beta2) * std::pow(R, 8) - 480 * std::pow(R, 4) * x(0) * ((-1.0 + beta2) * x(0) - 2 * x(1)) * std::pow(x(1), 2) + 1200 * (-1.0 + beta2) * std::pow(R, 2) * std::pow(x(0), 2) * std::pow(x(1), 4) + 50 * std::pow(R, 6) * x(1) * (-3 * x(0) + 2 * (-1.0 + beta2) * x(1)) + (47 * (-1.0 + beta2) * std::pow(R, 6) + 120 * std::pow(R, 2) * x(0) * (11 * (-1.0 + beta2) * x(0) - 27 * x(1)) * std::pow(x(1), 2) - 1680 * (-1.0 + beta2) * std::pow(x(0), 2) * std::pow(x(1), 4) + 10 * std::pow(R, 4) * x(1) * (81 * x(0) - 74 * (-1.0 + beta2) * x(1))) * r_sq + (-163 * (-1.0 + beta2) * std::pow(R, 4) - 840 * x(0) * ((-1.0 + beta2) * x(0) - 3 * x(1)) * std::pow(x(1), 2) + 10 * std::pow(R, 2) * x(1) * (-129 * x(0) + 142 * (-1.0 + beta2) * x(1))) * std::pow(r_sq, 2) + (197 * (-1.0 + beta2) * std::pow(R, 2) + 30 * x(1) * (21 * x(0) - 26 * (-1.0 + beta2) * x(1))) * std::pow(r_sq, 3) - 78 * (-1.0 + beta2) * std::pow(r_sq, 4)));

            return (2.0 * (R_sq - r_sq) / R_sq) * PolyMesh2D::Functional::RowVector((-1 + beta2) * R_sq * (3 + 2 * R_sq) * x(1) + (1 - beta2) * r_sq * (3 + 5 * std::pow(r_sq, 2) + 2 * R_sq * (5 + R_sq) - r_sq * (8 + 7 * R_sq)) * x(1) - 2 * x(0) * (-15 * std::pow(r_sq, 3) + R_sq * (3 + 2 * R_sq) + 2 * (-1 + beta2) * (3 + 6 * R_sq + std::pow(R_sq, 2)) * x(0) * x(1) + 5 * std::pow(r_sq, 2) * (4 + 3 * R_sq + 4 * (-1 + beta2) * x(0) * x(1)) - r_sq * (6 + 16 * R_sq + 3 * std::pow(R_sq, 2) + 8 * (-1 + beta2) * (3 + 2 * R_sq) * x(0) * x(1))), (-1 + beta2) * R_sq * (3 + 2 * R_sq) * x(0) + (1 - beta2) * r_sq * (3 + 5 * std::pow(r_sq, 2) + 2 * R_sq * (5 + R_sq) - r_sq * (8 + 7 * R_sq)) * x(0) - 2 * x(1) * (-15 * std::pow(r_sq, 3) + R_sq * (3 + 2 * R_sq) + 2 * (-1 + beta2) * (3 + 6 * R_sq + std::pow(R_sq, 2)) * x(0) * x(1) + 5 * std::pow(r_sq, 2) * (4 + 3 * R_sq + 4 * (-1 + beta2) * x(0) * x(1)) - r_sq * (6 + 16 * R_sq + 3 * std::pow(R_sq, 2) + 8 * (-1 + beta2) * (3 + 2 * R_sq) * x(0) * x(1))));
        }
        else
        {
            // r > R
            // return PolyMesh2D::Functional::RowVector(-6 * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * (-4 * std::pow(x(0), 3) + 5 * (-1 + beta1) * std::pow(x(0), 2) * x(1) - 4 * x(0) * std::pow(x(1), 2) + (-1 + beta1) * std::pow(x(1), 3) + std::pow(R, 2) * (2 * x(0) + x(1) - beta1 * x(1))), 6 * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * (-((-1 + beta1) * std::pow(x(0), 3)) + std::pow(R, 2) * ((-1 + beta1) * x(0) - 2 * x(1)) + 4 * std::pow(x(0), 2) * x(1) - 5 * (-1 + beta1) * x(0) * std::pow(x(1), 2) + 4 * std::pow(x(1), 3)));

            // return (r_sq - R_sq) * PolyMesh2D::Functional::RowVector(x(1) * (3.0 * (1.0 - beta1) * std::pow(R, 8) + 50 * std::pow(R, 6) * x(0) * (2 * (-1.0 + beta1) * x(0) - 3 * x(1)) + 1200 * (-1.0 + beta1) * std::pow(R, 2) * std::pow(x(0), 4) * std::pow(x(1), 2) + 480 * std::pow(R, 4) * std::pow(x(0), 2) * x(1) * (2 * x(0) + x(1) - beta1 * x(1)) + (47 * (-1.0 + beta1) * std::pow(R, 6) - 1680 * (-1.0 + beta1) * std::pow(x(0), 4) * std::pow(x(1), 2) + 10 * std::pow(R, 4) * x(0) * (-74 * (-1.0 + beta1) * x(0) + 81 * x(1)) - 120 * std::pow(R, 2) * std::pow(x(0), 2) * x(1) * (27 * x(0) - 11 * (-1.0 + beta1) * x(1))) * r_sq + (-163 * (-1.0 + beta1) * std::pow(R, 4) + 10 * std::pow(R, 2) * x(0) * (142 * (-1.0 + beta1) * x(0) - 129 * x(1)) + 840 * std::pow(x(0), 2) * x(1) * (3 * x(0) + x(1) - beta1 * x(1))) * std::pow(r_sq, 2) + (197 * (-1.0 + beta1) * std::pow(R, 2) + 30 * x(0) * (-26 * (-1.0 + beta1) * x(0) + 21 * x(1))) * std::pow(r_sq, 3) - 78 * (-1.0 + beta1) * std::pow(r_sq, 4)), x(0) * (-3 * (-1.0 + beta1) * std::pow(R, 8) - 480 * std::pow(R, 4) * x(0) * ((-1.0 + beta1) * x(0) - 2 * x(1)) * std::pow(x(1), 2) + 1200 * (-1.0 + beta1) * std::pow(R, 2) * std::pow(x(0), 2) * std::pow(x(1), 4) + 50 * std::pow(R, 6) * x(1) * (-3 * x(0) + 2 * (-1.0 + beta1) * x(1)) + (47 * (-1.0 + beta1) * std::pow(R, 6) + 120 * std::pow(R, 2) * x(0) * (11 * (-1.0 + beta1) * x(0) - 27 * x(1)) * std::pow(x(1), 2) - 1680 * (-1.0 + beta1) * std::pow(x(0), 2) * std::pow(x(1), 4) + 10 * std::pow(R, 4) * x(1) * (81 * x(0) - 74 * (-1.0 + beta1) * x(1))) * r_sq + (-163 * (-1.0 + beta1) * std::pow(R, 4) - 840 * x(0) * ((-1.0 + beta1) * x(0) - 3 * x(1)) * std::pow(x(1), 2) + 10 * std::pow(R, 2) * x(1) * (-129 * x(0) + 142 * (-1.0 + beta1) * x(1))) * std::pow(r_sq, 2) + (197 * (-1.0 + beta1) * std::pow(R, 2) + 30 * x(1) * (21 * x(0) - 26 * (-1.0 + beta1) * x(1))) * std::pow(r_sq, 3) - 78 * (-1.0 + beta1) * std::pow(r_sq, 4)));

            return (2.0 * (R_sq - r_sq) / R_sq) * PolyMesh2D::Functional::RowVector((-1 + beta1) * R_sq * (3 + 2 * R_sq) * x(1) + (1 - beta1) * r_sq * (3 + 5 * std::pow(r_sq, 2) + 2 * R_sq * (5 + R_sq) - r_sq * (8 + 7 * R_sq)) * x(1) - 2 * x(0) * (-15 * std::pow(r_sq, 3) + R_sq * (3 + 2 * R_sq) + 2 * (-1 + beta1) * (3 + 6 * R_sq + std::pow(R_sq, 2)) * x(0) * x(1) + 5 * std::pow(r_sq, 2) * (4 + 3 * R_sq + 4 * (-1 + beta1) * x(0) * x(1)) - r_sq * (6 + 16 * R_sq + 3 * std::pow(R_sq, 2) + 8 * (-1 + beta1) * (3 + 2 * R_sq) * x(0) * x(1))), (-1 + beta1) * R_sq * (3 + 2 * R_sq) * x(0) + (1 - beta1) * r_sq * (3 + 5 * std::pow(r_sq, 2) + 2 * R_sq * (5 + R_sq) - r_sq * (8 + 7 * R_sq)) * x(0) - 2 * x(1) * (-15 * std::pow(r_sq, 3) + R_sq * (3 + 2 * R_sq) + 2 * (-1 + beta1) * (3 + 6 * R_sq + std::pow(R_sq, 2)) * x(0) * x(1) + 5 * std::pow(r_sq, 2) * (4 + 3 * R_sq + 4 * (-1 + beta1) * x(0) * x(1)) - r_sq * (6 + 16 * R_sq + 3 * std::pow(R_sq, 2) + 8 * (-1 + beta1) * (3 + 2 * R_sq) * x(0) * x(1))));
        }
    };

    // std::function<double(Eigen::Vector2d)> u = [&alpha, &R](const Eigen::Vector2d &x) -> double
    // {
    //     if(std::abs(x.norm() - R) < 1E-14)
    //     {
    //         return 0.0;
    //     }
    //     if (x.norm() < R)
    //     {
    //         // r < R
    //         return std::pow(R * R - std::pow(x(0), 2) - std::pow(x(1), 2), 2) * (R * R * ((3.0 + alpha) * std::pow(x(0), 2) + (1.0 + 3.0 * alpha) * std::pow(x(1), 2)) - (std::pow(x(0), 2) + std::pow(x(1), 2)) * ((9.0 + alpha) * std::pow(x(0), 2) + (1.0 + 9.0 * alpha) * std::pow(x(1), 2)));
    //     }
    //     else
    //     {
    //         // r > R
    //         return std::pow(R * R - std::pow(x(0), 2) - std::pow(x(1), 2), 2) * (R * R * ((3.0 + alpha) * std::pow(x(1), 2) + (1.0 + 3.0 * alpha) * std::pow(x(0), 2)) - (std::pow(x(1), 2) + std::pow(x(0), 2)) * ((9.0 + alpha) * std::pow(x(1), 2) + (1.0 + 9.0 * alpha) * std::pow(x(0), 2)));
    //     }
    // };
    std::function<Eigen::Matrix2d(PolyMesh2D::CurvedMesh::Cell *)> diffusion = [&beta1, &beta2, &R](PolyMesh2D::CurvedMesh::Cell *cell) -> Eigen::Matrix2d
    {
        Eigen::Matrix2d K = Eigen::Matrix2d::Identity(2, 2);

        if (cell->center_mass().norm() < R)
        {
            // r < R
            // K(0, 0) = alpha;
            K(0, 1) = 1.0 - beta1;
            K(1, 0) = 1.0 - beta1;
        }
        else
        {
            // r > R
            // K(1, 1) = alpha;
            K(0, 1) = 1.0 - beta2;
            K(1, 0) = 1.0 - beta2;
        }

        return K;
    };
    // std::function<Eigen::Matrix2d(PolyMesh2D::CurvedMesh::Cell *)> diffusion = [&alpha, &R](PolyMesh2D::CurvedMesh::Cell *cell) -> Eigen::Matrix2d
    // {
    //     Eigen::Matrix2d K = Eigen::Matrix2d::Identity(2, 2);

    //     if (cell->center_mass().norm() < R)*
    //     {
    //         // r < R
    //         // K(0, 0) = alpha;
    //         K(0, 1) = 1.0 - 1E-3;
    //         K(1, 0) = 1.0 - 1E-3;
    //     }
    //     else
    //     {
    //         // r > R
    //         // K(1, 1) = alpha;
    //     }

    //     return K;
    // };

    // std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [&beta, &R](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    // {
    //     // return 1.0;
    //     return 24.0 * (std::pow(R, 4) + 2.0 * x.squaredNorm() * (3.0 * x.squaredNorm() + 4.0 * x(0) * x(1) * (1.0 - beta)) - 6.0 * std::pow(R, 2) * (x.squaredNorm() + x(0) * x(1) * (1.0 - beta)));
    // };

    // std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [&beta1, &beta2, &R](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [&a, &R](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    {
        // double r_sq = x.squaredNorm();
        // double R_sq = R * R;
        // if(std::abs(r_sq - R_sq) < 1E-14)
        // {
        //     return 0.0;
        // }
        // return (1.0 - r_sq) * std::pow(std::abs(R_sq - r_sq), 5.0 / 3.0);
        // if (std::abs(x.norm() - R) < 1E-12)
        // {
        //     if (cell->center_mass().norm() < R)
        //     {
        //         return 4.0 * std::pow(R_sq - r_sq, a - 2.0) * (std::pow(a, 2) * r_sq * (r_sq - 1.0) + std::pow(R_sq - r_sq, 2) + a * (2 * std::pow(r_sq, 2) + R_sq - 3 * r_sq * R_sq));
        //     }
        //     else
        //     {
        //         return 4.0 * std::pow(r_sq - R_sq, a - 2.0) * (std::pow(a, 2) * r_sq * (r_sq - 1.0) + std::pow(r_sq - R_sq, 2) + a * (2 * std::pow(r_sq, 2) + R_sq - 3 * r_sq * R_sq));
        //     }
        // }
        // if (x.norm() < R)
        // {
        //     return 4.0 * std::pow(R_sq - r_sq, a - 2.0) * (std::pow(a, 2) * r_sq * (r_sq - 1.0) + std::pow(R_sq - r_sq, 2) + a * (2 * std::pow(r_sq, 2) + R_sq - 3 * r_sq * R_sq));
        // }
        // else
        // {
        //     return 4.0 * std::pow(r_sq - R_sq, a - 2.0) * (std::pow(a, 2) * r_sq * (r_sq - 1.0) + std::pow(r_sq - R_sq, 2) + a * (2 * std::pow(r_sq, 2) + R_sq - 3 * r_sq * R_sq));
        // }
        // double r = x.norm();
        return 1.0;
        // return 12 * ((3 + beta1 * (-1 + beta2) - beta2) * std::pow(R, 4) + (17 - 5 * beta1 + 5 * (-1 + beta1) * beta2) * std::pow(x(0), 4) - 16 * (-2 + beta1 + beta2) * std::pow(x(0), 3) * x(1) + 6 * (7 - 3 * beta1 + 3 * (-1 + beta1) * beta2) * std::pow(x(0), 2) * std::pow(x(1), 2) - 16 * (-2 + beta1 + beta2) * x(0) * std::pow(x(1), 3) + (17 - 5 * beta1 + 5 * (-1 + beta1) * beta2) * std::pow(x(1), 4) + 6 * std::pow(R, 2) * ((-3 + beta1 + beta2 - beta1 * beta2) * std::pow(x(0), 2) + 2 * (-2 + beta1 + beta2) * x(0) * x(1) + (-3 + beta1 + beta2 - beta1 * beta2) * std::pow(x(1), 2)));

        // return 2.0 * (-78.0 * (-1.0 + beta1) * (-1.0 + beta2) * std::pow(r, 10) + 3.0 * (-1.0 + beta1) * (-1.0 + beta2) * std::pow(R, 10) - 2880.0 * (-1.0 + beta2) * std::pow(R, 2) * std::pow(x(0), 3) * std::pow(x(1), 3) * (std::pow(x(0), 2) - 2 * (-1.0 + beta1) * x(0) * x(1) + std::pow(x(1), 2)) + 25 * std::pow(R, 8) * ((-7 - 4 * beta1 * (-1.0 + beta2) + 4 * beta2) * std::pow(x(0), 2) + 12 * (-2 + beta1 + beta2) * x(0) * x(1) + (-7 - 4 * beta1 * (-1.0 + beta2) + 4 * beta2) * std::pow(x(1), 2)) + 600 * std::pow(R, 4) * std::pow(x(0), 2) * std::pow(x(1), 2) * ((-13 - 6 * beta1 * (-1.0 + beta2) + 6 * beta2) * std::pow(x(0), 2) + 14 * (-2 + beta1 + beta2) * x(0) * x(1) + (-13 - 6 * beta1 * (-1.0 + beta2) + 6 * beta2) * std::pow(x(1), 2)) - 120 * std::pow(R, 6) * x(0) * x(1) * ((-27 + 16 * beta1 + 11 * beta2) * std::pow(x(0), 2) + 2 * (-33 - 13 * beta1 * (-1.0 + beta2) + 13 * beta2) * x(0) * x(1) + (-27 + 16 * beta1 + 11 * beta2) * std::pow(x(1), 2)) + 5 * std::pow(r, 8) * (55 * (-1.0 + beta1) * (-1.0 + beta2) * std::pow(R, 2) - 3 * (73 - 52 * beta1 + 52 * (-1.0 + beta1) * beta2) * std::pow(x(0), 2) + 36 * (-20 + 7 * beta1 + 13 * beta2) * x(0) * x(1) - 3 * (73 - 52 * beta1 + 52 * (-1.0 + beta1) * beta2) * std::pow(x(1), 2)) + 10 * std::pow(r, 2) * (-5 * (-1.0 + beta1) * (-1.0 + beta2) * std::pow(R, 8) + 336 * (-1.0 + beta2) * std::pow(x(0), 3) * std::pow(x(1), 3) * (std::pow(x(0), 2) - 2 * (-1.0 + beta1) * x(0) * x(1) + std::pow(x(1), 2)) + 12 * std::pow(R, 4) * x(0) * x(1) * ((-121 + 70 * beta1 + 51 * beta2) * std::pow(x(0), 2) + (-292 + 117 * beta1 - 117 * (-1.0 + beta1) * beta2) * x(0) * x(1) + (-121 + 70 * beta1 + 51 * beta2) * std::pow(x(1), 2)) + 288 * std::pow(R, 2) * std::pow(x(0), 2) * std::pow(x(1), 2) * ((7 - 3 * beta1 + 3 * (-1.0 + beta1) * beta2) * std::pow(x(0), 2) + (15 - 8 * beta1 - 7 * beta2) * x(0) * x(1) + (7 - 3 * beta1 + 3 * (-1.0 + beta1) * beta2) * std::pow(x(1), 2)) + 12 * std::pow(R, 6) * ((11 - 7 * beta1 + 7 * (-1.0 + beta1) * beta2) * std::pow(x(0), 2) + (37 - 16 * beta1 - 21 * beta2) * x(0) * x(1) + (11 - 7 * beta1 + 7 * (-1.0 + beta1) * beta2) * std::pow(x(1), 2))) + 40 * std::pow(r, 6) * (-9 * (-1.0 + beta1) * (-1.0 + beta2) * std::pow(R, 4) + 3 * x(0) * x(1) * ((-75 + 42 * beta1 + 33 * beta2) * std::pow(x(0), 2) + (-178 + 73 * beta1 - 73 * (-1.0 + beta1) * beta2) * x(0) * x(1) + 3 * (-25 + 14 * beta1 + 11 * beta2) * std::pow(x(1), 2)) + std::pow(R, 2) * ((79 - 55 * beta1 + 55 * (-1.0 + beta1) * beta2) * std::pow(x(0), 2) - 3 * (-87 + 32 * beta1 + 55 * beta2) * x(0) * x(1) + (79 - 55 * beta1 + 55 * (-1.0 + beta1) * beta2) * std::pow(x(1), 2))) + 30 * std::pow(r, 4) * (7 * (-1.0 + beta1) * (-1.0 + beta2) * std::pow(R, 6) + 4 * std::pow(R, 2) * x(0) * x(1) * ((169 - 96 * beta1 - 73 * beta2) * std::pow(x(0), 2) + 4 * (101 - 41 * beta1 + 41 * (-1.0 + beta1) * beta2) * x(0) * x(1) + (169.0 - 96 * beta1 - 73 * beta2) * std::pow(x(1), 2)) + std::pow(R, 4) * ((-107 + 72 * beta1 - 72 * (-1.0 + beta1) * beta2) * std::pow(x(0), 2) + 4 * (-89 + 35 * beta1 + 54 * beta2) * x(0) * x(1) + (-107 + 72 * beta1 - 72 * (-1.0 + beta1) * beta2) * std::pow(x(1), 2)) + 28 * std::pow(x(0), 2) * std::pow(x(1), 2) * ((-15 + 6 * beta1 + 6 * beta2 - 6 * beta1 * beta2) * std::pow(x(0), 2) + 2 * (-16.0 + 9 * beta1 + 7 * beta2) * x(0) * x(1) + 3 * (-5 + 2 * beta1 + 2 * beta2 - 2 * beta1 * beta2) * std::pow(x(1), 2))));

        // return (-4 * (3 * (-17 - 5 * beta1 * (-1 + beta2) + 5 * beta2) * std::pow(x(0), 4) + 8 * (27 + 7 * beta1 * (-1 + beta2) - 7 * beta2) * std::pow(x(0), 6) + 15 * (-13 - 3 * beta1 * (-1 + beta2) + 3 * beta2) * std::pow(x(0), 8) + 48 * (-2 + beta1 + beta2) * std::pow(x(0), 3) * (1 - 5 * std::pow(x(0), 2) + 5 * std::pow(x(0), 4)) * x(1) - 6 * std::pow(x(0), 2) * (21 - 9 * beta2 + 20 * (-7 + 3 * beta2) * std::pow(x(0), 2) + 10 * (17 - 7 * beta2) * std::pow(x(0), 4) + beta1 * (-1 + beta2) * (9 - 60 * std::pow(x(0), 2) + 70 * std::pow(x(0), 4))) * std::pow(x(1), 2) + 48 * (-2 + beta1 + beta2) * (x(0) - 10 * std::pow(x(0), 3) + 15 * std::pow(x(0), 5)) * std::pow(x(1), 3) - 3 * (17 - 5 * beta2 + 40 * (-7 + 3 * beta2) * std::pow(x(0), 2) + 50 * (11 - 5 * beta2) * std::pow(x(0), 4) + 5 * beta1 * (-1 + beta2) * (1 - 24 * std::pow(x(0), 2) + 50 * std::pow(x(0), 4))) * std::pow(x(1), 4) + 240 * (-2 + beta1 + beta2) * x(0) * (-1 + 3 * std::pow(x(0), 2)) * std::pow(x(1), 5) + 4 * (2 * (27 + 7 * beta1 * (-1 + beta2) - 7 * beta2) - 15 * (17 + 7 * beta1 * (-1 + beta2) - 7 * beta2) * std::pow(x(0), 2)) * std::pow(x(1), 6) + 240 * (-2 + beta1 + beta2) * x(0) * std::pow(x(1), 7) + 15 * (-13 - 3 * beta1 * (-1 + beta2) + 3 * beta2) * std::pow(x(1), 8) + 2 * std::pow(R, 6) * ((3 + beta1 * (-1 + beta2) - beta2) * (-1 + 3 * std::pow(x(0), 2)) - 6 * (-2 + beta1 + beta2) * x(0) * x(1) + 3 * (3 + beta1 * (-1 + beta2) - beta2) * std::pow(x(1), 2)) + 3 * std::pow(R, 4) * (-3 + beta1 + beta2 - beta1 * beta2 + 36 * std::pow(x(0), 2) - 12 * beta1 * std::pow(x(0), 2) - 12 * beta2 * std::pow(x(0), 2) + 12 * beta1 * beta2 * std::pow(x(0), 2) - 51 * std::pow(x(0), 4) + 15 * beta1 * std::pow(x(0), 4) + 15 * beta2 * std::pow(x(0), 4) - 15 * beta1 * beta2 * std::pow(x(0), 4) + 24 * (-2 + beta1 + beta2) * x(0) * (-1 + 2 * std::pow(x(0), 2)) * x(1) - 6 * (2 * (-3 + beta1 + beta2 - beta1 * beta2) + 3 * (7 + 3 * beta1 * (-1 + beta2) - 3 * beta2) * std::pow(x(0), 2)) * std::pow(x(1), 2) + 48 * (-2 + beta1 + beta2) * x(0) * std::pow(x(1), 3) + 3 * (-17 - 5 * beta1 * (-1 + beta2) + 5 * beta2) * std::pow(x(1), 4)) + 6 * std::pow(R, 2) * (3 * (3 + beta1 * (-1 + beta2) - beta2) * std::pow(x(0), 2) + 3 * (-17 - 5 * beta1 * (-1 + beta2) + 5 * beta2) * std::pow(x(0), 4) + 2 * (27 + 7 * beta1 * (-1 + beta2) - 7 * beta2) * std::pow(x(0), 6) - 6 * (-2 + beta1 + beta2) * (x(0) - 8 * std::pow(x(0), 3) + 10 * std::pow(x(0), 5)) * x(1) + 3 * (3 - beta2 + 6 * (-7 + 3 * beta2) * std::pow(x(0), 2) + 10 * (7 - 3 * beta2) * std::pow(x(0), 4) + beta1 * (-1 + beta2) * (1 - 18 * std::pow(x(0), 2) + 30 * std::pow(x(0), 4))) * std::pow(x(1), 2) - 24 * (-2 + beta1 + beta2) * x(0) * (-2 + 5 * std::pow(x(0), 2)) * std::pow(x(1), 3) + 3 * (-17 + 5 * beta2 + 10 * (7 - 3 * beta2) * std::pow(x(0), 2) + 5 * beta1 * (-1 + beta2) * (-1 + 6 * std::pow(x(0), 2))) * std::pow(x(1), 4) - 60 * (-2 + beta1 + beta2) * x(0) * std::pow(x(1), 5) + 2 * (27 + 7 * beta1 * (-1 + beta2) - 7 * beta2) * std::pow(x(1), 6)))) / std::pow(R, 2);
    };

    // std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [&alpha, &R](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    // {
    //     return -2 * (1 + alpha * (6 + alpha)) * std::pow(R * R, 3) + 36 * (1 + alpha * (6 + alpha)) * std::pow(R * R, 2) * (std::pow(x(0), 2) + std::pow(x(1), 2)) -
    //            18 * R * R * ((5 + alpha * (38 + 5 * alpha)) * std::pow(x(0), 4) + 6 * (3 + alpha) * (1 + 3 * alpha) * std::pow(x(0), 2) * std::pow(x(1), 2) + (5 + alpha * (38 + 5 * alpha)) * std::pow(x(1), 4)) +
    //            8 * (std::pow(x(0), 2) + std::pow(x(1), 2)) * ((7 + alpha * (66 + 7 * alpha)) * std::pow(x(0), 4) + 2 * (19 + alpha * (42 + 19 * alpha)) * std::pow(x(0), 2) * std::pow(x(1), 2) + (7 + alpha * (66 + 7 * alpha)) * std::pow(x(1), 4));
    // };
    // std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    // {
    //     return x(0) * x(1);
    // };

    // ScalarFunction2D source(src);

    // std::function<double(Eigen::Vector2d)> u = [](const Eigen::Vector2d &x) -> double
    // {
    //     double r_hat = std::sqrt(std::pow(x(0) + 0.25, 2) + std::pow(x(1) + 0.25, 2));

    //     return std::pow(1.0 - r_hat, 2) * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1));
    // };

    // std::function<double(Eigen::Vector2d)> u = [&alpha1, &alpha2, &level_set, &R](const Eigen::Vector2d &x) -> double
    // {
    //     // if(std::abs(level_set.value(x)) < 1E-12)
    //     // {
    //     //     return 0.0;
    //     // }
    //     double alpha = -1.0;
    //     // if(level_set.value(x) > 0.0)
    //     if(x.norm() < R)
    //     {
    //         // r < R
    //         alpha = alpha1;
    //     }
    //     else
    //     {
    //         // r > R
    //         alpha = alpha2;
    //     }
    //     // return level_set.value(x) * x(0) * x(1) / (6.0 * (1.0 + alpha));
    //     return (R * R - x.squaredNorm()) * x(0) * x(1) / (6.0 * (1.0 + alpha));
    // };

    // std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Du = [](const Eigen::Vector2d &x) -> PolyMesh2D::Functional::RowVector
    // {
    //     double r_hat = std::sqrt(std::pow(x(0) + 0.25, 2) + std::pow(x(1) + 0.25, 2));

    //     return (1.0 - r_hat) * PolyMesh2D::Functional::RowVector(std::sin(Math::PI * x(1)) * (Math::PI * (1.0 - r_hat) * std::cos(Math::PI * x(0)) - 2.0 * (0.25 + x(0)) * std::sin(Math::PI * x(0)) / r_hat), std::sin(Math::PI * x(0)) * (Math::PI * (1.0 - r_hat) * std::cos(Math::PI * x(1)) - 2.0 * (0.25 + x(1R)) * std::sin(Math::PI * x(1)) / r_hat));
    // };

    // std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Du = [](const Eigen::Vector2d &x) -> PolyMesh2D::Functional::RowVector
    // {
    //     return PolyMesh2D::Functional::RowVector(std::pow(-1 + std::pow(0.25 + x(0), 2) + std::pow(0.25 + x(1), 2), 2) * (Math::PI * (1 - std::pow(0.25 + x(0), 2) - std::pow(0.25 + x(1), 2)) * std::cos(Math::PI * x(0)) - 6 * (0.25 + x(0)) * std::sin(Math::PI * x(0))) * std::sin(Math::PI * x(1)), std::pow(-1 + std::pow(0.25 + x(0), 2) + std::pow(0.25 + x(1), 2), 2) * std::sin(Math::PI * x(0)) * (Math::PI * (1 - std::pow(0.25 + x(0), 2) - std::pow(0.25 + x(1), 2)) * std::cos(Math::PI * x(1)) - 6 * (0.25 + x(1)) * std::sin(Math::PI * x(1))));
    // };

    // std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Du = [](const Eigen::Vector2d &x) -> PolyMesh2D::Functional::RowVector
    // {
    //     return PolyMesh2D::Functional::RowVector(1.0, 1.0);
    // };

    // std::function<PolyMesh2D::Functional::RowVector(Eigen::Vector2d)> Du = [&alpha, &R](const Eigen::Vector2d &x) -> PolyMesh2D::Functional::RowVector
    // {
    //     if(std::abs(x.norm() - R) < 1E-14)
    //     {
    //         return PolyMesh2D::Functional::RowVector(0.0, 0.0);
    //     }
    //     if (x.norm() < R)
    //     {
    //         // r < R
    //         return PolyMesh2D::Functional::RowVector(2.0 * x(0) * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * ((3.0 + alpha) * std::pow(R, 4) + 4.0 * (std::pow(x(0), 2) + std::pow(x(1), 2)) * ((9.0 + alpha) * std::pow(x(0), 2) + (3.0 + 7.0 * alpha) * std::pow(x(1), 2)) - std::pow(R, 2) * ((27.0 + 5.0 * alpha) * std::pow(x(0), 2) + (15.0 + 17.0 * alpha) * std::pow(x(1), 2))), 2.0 * x(1) * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * ((1.0 + 3.0 * alpha) * std::pow(R, 4) + 4.0 * (std::pow(x(0), 2) + std::pow(x(1), 2)) * ((7.0 + 3.0 * alpha) * std::pow(x(0), 2) + (1.0 + 9.0 * alpha) * std::pow(x(1), 2)) - std::pow(R, 2) * ((17.0 + 15.0 * alpha) * std::pow(x(0), 2) + (5.0 + 27.0 * alpha) * std::pow(x(1), 2))));
    //     }
    //     else
    //     {
    //         // r > R
    //         return PolyMesh2D::Functional::RowVector(2.0 * x(0) * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * ((1.0 + 3.0 * alpha) * std::pow(R, 4) + 4.0 * (std::pow(x(0), 2) + std::pow(x(1), 2)) * ((1.0 + 9.0 * alpha) * std::pow(x(0), 2) + (7.0 + 3.0 * alpha) * std::pow(x(1), 2)) - std::pow(R, 2) * ((5.0 + 27.0 * alpha) * std::pow(x(0), 2) + (17.0 + 15.0 * alpha) * std::pow(x(1), 2))), 2.0 * x(1) * (std::pow(R, 2) - std::pow(x(0), 2) - std::pow(x(1), 2)) * ((3.0 + alpha) * std::pow(R, 4) + 4.0 * (std::pow(x(0), 2) + std::pow(x(1), 2)) * ((3.0 + 7.0 * alpha) * std::pow(x(0), 2) + (9.0 + alpha) * std::pow(x(1), 2)) - std::pow(R, 2) * ((15.0 + 17.0 * alpha) * std::pow(x(0), 2) + (27.0 + 5.0 * alpha) * std::pow(x(1), 2))));
    //     }
    // };

    ScalarFunction2D exact_sol(u, Du);

    // std::function<Eigen::Matrix2d(PolyMesh2D::CurvedMesh::Cell *)> diffusion = [&anisotropy, &level_set](PolyMesh2D::CurvedMesh::Cell *cell) -> Eigen::Matrix2d
    // {
    // std::function<Eigen::Matrix2d(PolyMesh2D::CurvedMesh::Cell *)> diffusion = [&alpha1, &alpha2, &level_set, &R](PolyMesh2D::CurvedMesh::Cell *cell) -> Eigen::Matrix2d
    // {
    //     Eigen::Matrix2d K = Eigen::Matrix2d::Identity(2, 2);
    //     // if (level_set.value(cell->center_mass()) > 0.0)
    //     // {
    //     //     // std::cout << "hello\n";
    //     //     K(0, 0) = anisotropy;
    //     // }
    //     // else
    //     // {
    //     //     std::cout << "hello2\n";
    //     // }

    //     double alpha = -1.0;
    //     // if(level_set.value(cell->center_mass()) > 0.0)
    //     if(cell->center_mass().norm() < R)
    //     {
    //         // r < R
    //         alpha = alpha1;
    //     }
    //     else
    //     {
    //         // r > R
    //         alpha = alpha2;
    //     }

    //     K(0, 0) = alpha;

    //     return K;
    // };

    // Model model(hho, source);
    Model model(hho, src, diffusion);

    size_t n_fixed_dofs = 0;
    for (auto &b_edge : curved_mesh->get_b_edges())
    {
        n_fixed_dofs += hho.local_edge_dofs(b_edge->global_index());
    }

    // std::cout << curved_mesh->n_cells() << "\n";

    // double integral = 0.0;
    // for (size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
    // {
    //     Quadrature::QuadratureRule<Eigen::Vector2d> cell_quad = hho.quadT(iT);
    //     for(size_t iqn = 0; iqn < cell_quad.size(); ++iqn)
    //     {
    //         integral += cell_quad[iqn].w;
    //     }
    // }

    // std::cout << integral << "\n";
    // exit(0);

    // Eigen::VectorXd UDir = hho.interpolate(u).tail(n_fixed_dofs);

    // for (auto &b_edge : curved_mesh->get_b_edges())
    // {
    //     size_t iE = b_edge->global_index();
    //     size_t offset = hho.global_offset_E(iE) - hho.global_offset_E(curved_mesh->n_i_edges());
    //     Quadrature::QuadratureRule<double> edge_quad = hho.quadE(iE);

    //     for (size_t iqn = 0; iqn < edge_quad.size(); ++iqn)
    //     {
    //         for (size_t i = 0; i < hho.local_edge_dofs(iE); ++i)
    //         {
    //             UDir(offset + i) += edge_quad[iqn].w * hho.edge_basis(iE)->value(i, edge_quad[iqn].x) * u(b_edge->parameterisation().value(edge_quad[iqn].x));
    //         }
    //     }
    // }

    model.assemble(params.use_threads);

    Eigen::VectorXd UDir(Eigen::VectorXd::Zero(n_fixed_dofs));
    Eigen::VectorXd approx_Uvec(model.solve(UDir));

    // std::cout << "\n[Scheme] Computing interpolant\n";
    // Eigen::VectorXd interp_Uvec(hho.interpolate(u));
    size_t n_dofs = hho.total_cell_dofs() + hho.total_edge_dofs();
    Eigen::VectorXd interp_Uvec(Eigen::VectorXd::Zero(n_dofs));

    auto errors = model.compute_errors(approx_Uvec, interp_Uvec, exact_sol);
    // double l2_error = errors[0];
    // double h1_error = errors[1];
    // double energy_error = errors[2];

    // double l2_error = errors[0];
    // double h1_error = errors[1];
    // double energy_error = errors[0];

    double integral = errors[0];
    double integral_error = errors[1];
    double h1_norm = errors[2];
    double h1_norm_error = errors[3];

    // std::cout << "     L2 Error = " << l2_error << "\n";
    // std::cout << "     H1 Error = " << h1_error << "\n";
    // std::cout << "     Energy Error = " << energy_error << "\n";

    model.plot(approx_Uvec, interp_Uvec, params.plot_file, exact_sol);

    std::cout << "\n[Scheme] Writing solution to file\n";

    std::ofstream out("results.txt");
    out << "EdgeDegree: " << params.edge_degree << "\n";
    out << "CellDegree: " << params.cell_degree << "\n";
    // out << std::setprecision(20) << "L2Error: " << l2_error << "\n";
    // out << std::setprecision(20) << "H1Error: " << h1_error << "\n";
    // out << std::setprecision(20) << "EnergyError: " << energy_error << "\n";
    out << std::setprecision(20) << "integral: " << integral << "\n";
    out << std::setprecision(20) << "integral_error: " << integral_error << "\n";
    out << std::setprecision(20) << "h1_norm: " << h1_norm << "\n";
    out << std::setprecision(20) << "h1_norm_error: " << h1_norm_error << "\n";
    out << "MeshSize: " << curved_mesh->h_max() << "\n";
    out << "NbCells: " << curved_mesh->n_cells() << "\n";
    out << "NbEdges: " << curved_mesh->n_edges() << "\n";
    out << "NbInternalEdges: " << curved_mesh->n_edges() - curved_mesh->n_b_edges() << "\n";
    out.close();

    return 0;
}

ModelParameters::ModelParameters(const int argc, const char **argv)
{
    namespace po = boost::program_options;

    // Program options
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Produce help message")("mesh,m", po::value<std::string>(), "Set the mesh")("celldegree,l", po::value<unsigned>(), "Set the degree of the cell polynomials")("edgedegree,k", po::value<unsigned>(), "Set the degree of the edge polynomials")("plot,p", po::value<std::string>(), "Plot to file")("use_threads,u", po::value<bool>(), "Using multithreading")("orthonormalise,o", po::value<bool>(), "Orthonormalise the basis functions");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Display help message
    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        exit(0);
    }

    // Get mesh file
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "/home/liam/github/codes/liamyemm/HArDCore/2D/typ2_meshes/mesh2_2.typ2");

    // Get polynomial degrees
    edge_degree = (vm.count("edgedegree") ? vm["edgedegree"].as<unsigned>() : 0);
    cell_degree = (vm.count("celldegree") ? vm["celldegree"].as<unsigned>() : edge_degree);

    if ((std::abs((int)cell_degree - (int)edge_degree) > 1))
    {
        std::cout << "Invalid choice of polynomial degrees\n";
        exit(1);
    }

    // Get plot file
    plot_file = (vm.count("plot") ? vm["plot"].as<std::string>() : "");

    // Get use_threads
    use_threads = (vm.count("use_threads") ? vm["use_threads"].as<bool>() : true);

    orthonormalise = (vm.count("orthonormalise") ? vm["orthonormalise"].as<bool>() : true);
}

Model::Model(const HybridCore &hho, const std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> &src, std::function<Eigen::Matrix2d(CurvedMesh::Cell *)> diffusion) : m_hho(hho), m_src(src), m_diffusion(diffusion), mesh_ptr(m_hho.get_mesh()) {}

void Model::assemble(bool use_threads)
{
    std::cout << "\n[Scheme] Assembling the linear system\n";

    // Start the timer
    boost::timer::cpu_timer assembly_timer;
    assembly_timer.start();

    // Set up triplets for sparse matrix initialisation
    std::vector<Eigen::Triplet<double>> triplets_GlobMat;
    std::vector<Eigen::Triplet<double>> triplets_ScBe;

    // set up vectors of local matrices to emplace into global matrices
    std::vector<Eigen::MatrixXd> invATT_ATF(mesh_ptr->n_cells());
    std::vector<Eigen::VectorXd> cell_source(mesh_ptr->n_cells());
    std::vector<Eigen::MatrixXd> MatF(mesh_ptr->n_cells());

    GlobRHS = Eigen::VectorXd::Zero(m_hho.total_edge_dofs());
    ScRHS = Eigen::VectorXd::Zero(m_hho.total_cell_dofs());

    GlobMat.resize(m_hho.total_edge_dofs(), m_hho.total_edge_dofs());
    ScBeMat.resize(m_hho.total_cell_dofs(), m_hho.total_edge_dofs());

    std::vector<Eigen::VectorXd> bT;

    PT.resize(mesh_ptr->n_cells());
    AT.resize(mesh_ptr->n_cells());
    bT.resize(mesh_ptr->n_cells());

    // Construct the local matrices using multithreading if use_threads is true
    std::function<void(size_t, size_t)> construct_all_local_contributions = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        // for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
        {
            const size_t n_local_bdry_dofs = m_hho.local_boundary_dofs(iT);
            const size_t n_local_cell_dofs = m_hho.local_cell_dofs(iT);

            local_poisson_operator(iT);
            local_source_term(iT, bT[iT]);

            // Get each contribution of AT
            Eigen::MatrixXd ATT(AT[iT].topLeftCorner(n_local_cell_dofs, n_local_cell_dofs));
            Eigen::MatrixXd ATF(AT[iT].topRightCorner(n_local_cell_dofs, n_local_bdry_dofs));
            Eigen::MatrixXd AFT(AT[iT].bottomLeftCorner(n_local_bdry_dofs, n_local_cell_dofs));
            Eigen::MatrixXd AFF(AT[iT].bottomRightCorner(n_local_bdry_dofs, n_local_bdry_dofs));

            Eigen::PartialPivLU<Eigen::MatrixXd> invATT;
            invATT.compute(ATT);

            Eigen::VectorXd invATT_bTcell = invATT.solve(bT[iT].head(n_local_cell_dofs));

            // Store the local matrices
            invATT_ATF[iT] = invATT.solve(ATF);
            MatF[iT] = AFF - AFT * invATT_ATF[iT];
            cell_source[iT] = bT[iT].tail(n_local_bdry_dofs) - AFT * invATT_bTcell;

            ScRHS.segment(m_hho.global_offset_T(iT), n_local_cell_dofs) = invATT_bTcell;
        }
    };

    // Running the local constructions in parallel
    parallel_for(mesh_ptr->n_cells(), construct_all_local_contributions, use_threads);

    // Set the local matrices into the correct positions of the global matrices
    for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++)
    {
        const size_t n_local_cell_dofs = m_hho.local_cell_dofs(iT);
        PolyMesh2D::CurvedMesh::Cell *cell = mesh_ptr->cell(iT);
        for (size_t iTF = 0; iTF < cell->n_edges(); iTF++)
        {
            const size_t iF = cell->edge(iTF)->global_index();
            const size_t global_offset_iF = m_hho.global_offset_E(iF);
            const size_t local_offset_iTF = m_hho.local_offset_E(iT, iTF);
            for (size_t ik = 0; ik < m_hho.local_edge_dofs(iF); ik++)
            {
                const size_t iLocal = local_offset_iTF + ik;
                const size_t iGlobal = global_offset_iF + ik;
                GlobRHS(iGlobal) += cell_source[iT](iLocal);
                for (size_t i = 0; i < n_local_cell_dofs; i++)
                {
                    triplets_ScBe.emplace_back(m_hho.global_offset_T(iT) + i, iGlobal, invATT_ATF[iT](i, iLocal));
                }
                for (size_t jTF = 0; jTF < cell->n_edges(); jTF++)
                {
                    const size_t jF = cell->edge(jTF)->global_index();
                    const size_t global_offset_jF = m_hho.global_offset_E(jF);
                    const size_t local_offset_jTF = m_hho.local_offset_E(iT, jTF);
                    for (size_t jk = 0; jk < m_hho.local_edge_dofs(jF); jk++)
                    {
                        const size_t jLocal = local_offset_jTF + jk;
                        const size_t jGlobal = global_offset_jF + jk;
                        triplets_GlobMat.emplace_back(iGlobal, jGlobal, MatF[iT](iLocal, jLocal));
                    }
                }
            }
        }
    }

    // Construct the global matrices
    GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));
    ScBeMat.setFromTriplets(std::begin(triplets_ScBe), std::end(triplets_ScBe));

    std::cout << "     Assembly time = " << assembly_timer.elapsed().wall * pow(10, -9) << "s\n";
}

Eigen::VectorXd Model::solve(const Eigen::VectorXd &UDir)
{
    std::cout << "\n[Scheme] Solving the linear system\n";

    // Start the timer
    boost::timer::cpu_timer timer;
    timer.start();

    size_t n_fixed_dofs = UDir.size();

    // for (size_t iF = 0; iF < mesh_ptr->n_edges(); ++iF)
    // {
    //     if (!(mesh_ptr->edge(iF)->is_boundary()))
    //     {
    //         continue;
    //     }
    //     n_fixed_dofs += m_hho.local_edge_dofs(iF);
    // }
    const size_t n_unknowns = m_hho.total_edge_dofs() - n_fixed_dofs;

    Eigen::SparseMatrix<double> SysMat(GlobMat.topLeftCorner(n_unknowns, n_unknowns));
    // Eigen::VectorXd UDir = Eigen::VectorXd::Zero(n_fixed_dofs);

    Eigen::VectorXd B(GlobRHS.segment(0, n_unknowns) - GlobMat.topRightCorner(n_unknowns, n_fixed_dofs) * UDir);

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    // Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(SysMat);

    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    // solver.analyzePattern(SysMat);
    // solver.factorize(SysMat);

    // xF = INV( AII - AIT * INV_ATT * ATI ) * ( BI - AIT * INV_ATT * BT - (AID - AIT * INV_ATT * ATD * UD) )
    Eigen::VectorXd xF(solver.solve(B));

    // Print solver iterations and estimated error
    std::cout << "     [solver] #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;
    // _iterations = solver.iterations();

    // Find the residual of the system
    std::cout << "     [solver] residual: " << (SysMat * xF - B).norm() / B.norm() << std::endl;

    // Set each part of the solution
    Eigen::VectorXd Xh(Eigen::VectorXd::Zero(m_hho.total_cell_dofs() + m_hho.total_edge_dofs()));
    Xh.tail(n_fixed_dofs) = UDir;
    Xh.segment(m_hho.total_cell_dofs(), n_unknowns) = xF;

    //                         = INV_ATT_BT - INV_ATT * ATF * UF
    Xh.head(m_hho.total_cell_dofs()) = ScRHS - ScBeMat * Xh.tail(m_hho.total_edge_dofs());

    std::cout << "     Solving time = " << timer.elapsed().wall * pow(10, -9) << "s\n";

    return Xh;
}

std::vector<double> Model::compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const ScalarFunction2D &sol)
{
    std::vector<double> pT_average_vec(mesh_ptr->n_cells());

    std::vector<double> pT_H1_norm_vec(mesh_ptr->n_cells());
    // std::vector<double> flux_vec(mesh_ptr->n_cells());
    // double flux = 0.0;

    std::function<void(size_t, size_t)> compute_average = [&](size_t start, size_t end) -> void
    {
    for (size_t iT = start; iT < end; iT++)
    {
    // for (size_t iiF = 0; iiF < mesh_ptr->n_i_edges(); ++iiF)
    // {
        // size_t iF = mesh_ptr->i_edge(iiF)->global_index();
        // size_t iT1 = mesh_ptr->i_edge(iiF)->cell(0)->global_index();
        // size_t iT2 = mesh_ptr->i_edge(iiF)->cell(1)->global_index();

        // double HT = mesh_ptr->i_edge(iiF)->cell(0)->diam() + mesh_ptr->i_edge(iiF)->cell(1)->diam();
        // // QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

        // Eigen::VectorXd pT1_vec = PT[iT1] * m_hho.restr(approx_Uvec, iT1);
        // auto basis1 = m_hho.highorder_basis(iT1);

        // Eigen::VectorXd pT2_vec = PT[iT2] * m_hho.restr(approx_Uvec, iT2);
        // auto basis2 = m_hho.highorder_basis(iT2);

        Eigen::VectorXd pT_vec = PT[iT] * m_hho.restr(approx_Uvec, iT);
        auto basis = m_hho.highorder_basis(iT);

        std::function<double(PolyMesh2D::Functional::ColVector)> approx_sol = [&pT_vec, &basis](const PolyMesh2D::Functional::ColVector &x) -> double
        {
            double value = 0.0;
            for (size_t i = 0; i < basis->dimension(); i++)
            {
                value += pT_vec(i) * basis->value(i, x);
            }
            return value;
        };

        std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> approx_grad = [&pT_vec, &basis](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
        {
            PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
            for (size_t i = 0; i < basis->dimension(); i++)
            {
                value += pT_vec(i) * basis->derivative(i, x);
            }
            return value;
        };

        pT_average_vec[iT] = 0.0;
        pT_H1_norm_vec[iT] = 0.0;
        QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

        for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
        {
            pT_average_vec[iT] += quadT[iqn].w * approx_sol(quadT[iqn].x);
            pT_H1_norm_vec[iT] += quadT[iqn].w * approx_grad(quadT[iqn].x) * approx_grad(quadT[iqn].x).transpose();
        }

        // Eigen::Matrix2d diff_on_cell = m_diffusion(mesh_ptr->cell(iT));
        // Eigen::Matrix2d diff_T1 = m_diffusion(mesh_ptr->cell(iT1));
        // Eigen::Matrix2d diff_T2 = m_diffusion(mesh_ptr->cell(iT2));
        // std::function<Eigen::Matrix2d(Eigen::Vector2d)> diff_T1 = [&](Eigen::Vector2d x) -> Eigen::Matrix2d
        // {
        //     Eigen::Matrix2d K = Eigen::Matrix2d::Identity(2, 2);
        //     double beta1 = 1E-6;
        //     double beta2 = 1.0;
        //     double R = 0.8;
        //     if (std::abs(x.norm() - R) < 1E-14)
        //     {
        //         if (mesh_ptr->cell(iT1)->center_mass().norm() < R)
        //         {
        //             // r < R
        //             K(0, 1) = 1.0 - beta1;
        //             K(1, 0) = 1.0 - beta1;
        //         }
        //         else
        //         {
        //             // r > R
        //             K(0, 1) = 1.0 - beta2;
        //             K(1, 0) = 1.0 - beta2;
        //         }
        //     }
        //     else
        //     {
        //         if (x.norm() < R)
        //         {
        //             // r < R
        //             K(0, 1) = 1.0 - beta1;
        //             K(1, 0) = 1.0 - beta1;
        //         }
        //         else
        //         {
        //             // r > R
        //             K(0, 1) = 1.0 - beta2;
        //             K(1, 0) = 1.0 - beta2;
        //         }
        //     }

        //     return K;
        // };
        // std::function<Eigen::Matrix2d(Eigen::Vector2d)> diff_T2 = [&](Eigen::Vector2d x) -> Eigen::Matrix2d
        // {
        //     Eigen::Matrix2d K = Eigen::Matrix2d::Identity(2, 2);
        //     double beta1 = 1E-6;
        //     double beta2 = 1.0;
        //     double R = 0.8;
        //     if (std::abs(x.norm() - R) < 1E-12)
        //     {
        //         if (mesh_ptr->cell(iT2)->center_mass().norm() < R)
        //         {
        //             // r < R
        //             K(0, 1) = 1.0 - beta1;
        //             K(1, 0) = 1.0 - beta1;
        //         }
        //         else
        //         {
        //             // r > R
        //             K(0, 1) = 1.0 - beta2;
        //             K(1, 0) = 1.0 - beta2;
        //         }
        //     }
        //     else
        //     {
        //         if (x.norm() < R)
        //         {
        //             // r < R
        //             K(0, 1) = 1.0 - beta1;
        //             K(1, 0) = 1.0 - beta1;
        //         }
        //         else
        //         {
        //             // r > R
        //             K(0, 1) = 1.0 - beta2;
        //             K(1, 0) = 1.0 - beta2;
        //         }
        //     }

        //     return K;
        // };

        // std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> diff_grad_pT1 = [&pT1_vec, &basis1, &diff_T1](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
        // {
        //     PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
        //     for (size_t i = 0; i < basis1->dimension(); i++)
        //     {
        //         value += pT1_vec(i) * basis1->derivative(i, x) * diff_T1(x);
        //     }
        //     return value;
        // };

        // std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> diff_grad_pT2 = [&pT2_vec, &basis2, &diff_T2](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
        // {
        //     PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
        //     for (size_t i = 0; i < basis2->dimension(); i++)
        //     {
        //         value += pT2_vec(i) * basis2->derivative(i, x) * diff_T2(x);
        //     }
        //     return value;
        // };

        // std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> diff_grad_pT1 = [&pT1_vec, &basis1, &diff_exact](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
        // {
        //     PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
        //     for (size_t i = 0; i < basis1->dimension(); i++)
        //     {
        //         value += pT1_vec(i) * basis1->derivative(i, x) * diff_exact(x);
        //     }
        //     return value;
        // };

        // std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> diff_grad_pT2 = [&pT2_vec, &basis2, &diff_exact](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
        // {
        //     PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
        //     for (size_t i = 0; i < basis2->dimension(); i++)
        //     {
        //         value += pT2_vec(i) * basis2->derivative(i, x) * diff_exact(x);
        //     }
        //     return value;
        // };

        // flux_vec[iT] = 0.0;
        // for(size_t iTF = 0; iTF < mesh_ptr->cell(iT)->n_edges(); ++iTF)
        // {
        //     if(mesh_ptr->cell(iT)->edge(iTF)->is_boundary())
        //     {
        //         continue;
        //     }
        // QuadratureRule<double> quadF = m_hho.quadE(iF);
        // PolyMesh2D::Functional::Curve edge_param(mesh_ptr->edge(iF)->parameterisation());

        // double face_flux = 0.0;

        // for (size_t iqnF = 0; iqnF < quadF.size(); ++iqnF)
        // {
        //     // face_flux += quadF[iqnF].w * diff_grad_pT1(edge_param.value(quadF[iqnF].x)) * mesh_ptr->cell(iT)->edge_normal(iTF, quadF[iqnF].x);
        //     face_flux += quadF[iqnF].w * std::pow((diff_grad_pT1(edge_param.value(quadF[iqnF].x)) - diff_grad_pT2(edge_param.value(quadF[iqnF].x))) * edge_param.normal(quadF[iqnF].x), 2);
        // }

        // flux += HT * face_flux;

        // }
    }
    };
    parallel_for(mesh_ptr->n_cells(), compute_average, true);
    double pT_average = 0.0;
    double pT_H1_norm = 0.0;
    // double flux = 0.0;

    for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
    {
        pT_average += pT_average_vec[iT];
        pT_H1_norm += pT_H1_norm_vec[iT];
        // flux += flux_vec[iT];
    }
    // return {std::abs(pT_average - 0.47342316790031907514), std::abs(flux + Math::PI)};
    // return {pT_average, std::abs(pT_average - 0.46006947132345432649), std::sqrt(pT_H1_norm), std::abs(std::sqrt(pT_H1_norm) - 0.80699765696519920599)};
    return {pT_average, std::abs(pT_average - 0.46006947132349923502), std::sqrt(pT_H1_norm), std::abs(std::sqrt(pT_H1_norm) - 0.8069976569614020212)};
    
    // return {pT_average, std::sqrt(pT_H1_norm)};
    // return {std::sqrt(flux)};

    std::cout << "\n[Scheme] Computing error\n";

    std::vector<double> energy_diff_vec(mesh_ptr->n_cells());
    std::vector<double> energy_exact_vec(mesh_ptr->n_cells());

    std::vector<double> L2_diff_vec(mesh_ptr->n_cells());
    std::vector<double> L2_exact_vec(mesh_ptr->n_cells());

    std::vector<double> H1_diff_vec(mesh_ptr->n_cells());
    std::vector<double> H1_exact_vec(mesh_ptr->n_cells());

    // double H1_exact = 1.0;

    std::function<void(size_t, size_t)> compute_all_errors = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            Eigen::VectorXd ITk = m_hho.restr(interp_Uvec, iT);
            Eigen::VectorXd diff_vec = m_hho.restr(approx_Uvec, iT) - ITk;

            energy_diff_vec[iT] = diff_vec.transpose() * AT[iT] * diff_vec;
            energy_exact_vec[iT] = ITk.transpose() * AT[iT] * ITk;

            QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

            Eigen::VectorXd pT_vec = PT[iT] * m_hho.restr(approx_Uvec, iT);
            auto basis = m_hho.highorder_basis(iT);

            Eigen::Matrix2d diff_on_cell = m_diffusion(mesh_ptr->cell(iT));

            std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> approx_grad = [&pT_vec, &basis](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
            {
                PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
                for (size_t i = 0; i < basis->dimension(); i++)
                {
                    value += pT_vec(i) * basis->derivative(i, x);
                }
                return value;
            };

            H1_diff_vec[iT] = 0.0;
            H1_exact_vec[iT] = 0.0;

            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                H1_diff_vec[iT] += quadT[iqn].w * ((approx_grad(quadT[iqn].x) - sol.derivative(quadT[iqn].x)) * diff_on_cell) * ((approx_grad(quadT[iqn].x) - sol.derivative(quadT[iqn].x)).transpose());
                H1_exact_vec[iT] += quadT[iqn].w * (sol.derivative(quadT[iqn].x) * diff_on_cell) * (sol.derivative(quadT[iqn].x).transpose());
            }

            std::function<double(PolyMesh2D::Functional::ColVector)> approx_sol = [&pT_vec, &basis](const PolyMesh2D::Functional::ColVector &x) -> double
            {
                double value = 0.0;
                for (size_t i = 0; i < basis->dimension(); i++)
                {
                    value += pT_vec(i) * basis->value(i, x);
                }
                return value;
            };

            L2_diff_vec[iT] = 0.0;
            L2_exact_vec[iT] = 0.0;

            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                L2_diff_vec[iT] += quadT[iqn].w * std::pow(approx_sol(quadT[iqn].x) - sol.value(quadT[iqn].x), 2);
                L2_exact_vec[iT] += quadT[iqn].w * std::pow(sol.value(quadT[iqn].x), 2);
            }
        }
    };
    parallel_for(mesh_ptr->n_cells(), compute_all_errors, true);

    double energy_diff = 0.0;
    double energy_exact = 0.0;

    double L2_diff = 0.0;
    double L2_exact = 0.0;

    double H1_diff = 0.0;
    double H1_exact = 0.0;

    for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
    {
        energy_diff += energy_diff_vec[iT];
        energy_exact += energy_exact_vec[iT];
        L2_diff += L2_diff_vec[iT];
        L2_exact += L2_exact_vec[iT];
        H1_diff += H1_diff_vec[iT];
        H1_exact += H1_exact_vec[iT];
    }

    return {std::sqrt(L2_diff / L2_exact), std::sqrt(H1_diff / H1_exact), std::sqrt(energy_diff / energy_exact)};
}

void Model::plot(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::string &plot_file, const ScalarFunction2D &sol)
{
    std::ofstream sol_out("solution_plot.tsv");

    for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++)
    {
        QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

        Eigen::VectorXd pT_vec = PT[iT] * m_hho.restr(approx_Uvec, iT);
        auto basis = m_hho.highorder_basis(iT);

        std::function<double(PolyMesh2D::Functional::ColVector)> approx_sol = [&pT_vec, &basis](const PolyMesh2D::Functional::ColVector &x) -> double
        {
            double value = 0.0;
            for (size_t i = 0; i < basis->dimension(); i++)
            {
                value += pT_vec(i) * basis->value(i, x);
            }
            return value;
        };

        for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
        {
            Eigen::Vector2d point(quadT[iqn].x);
            sol_out << std::setprecision(12) << std::setw(20) << std::left << point(0) << std::setprecision(12) << std::setw(20) << std::left << point(1) << std::setprecision(12) << std::setw(20) << approx_sol(point) << std::endl;
        }
    }
    sol_out.close();

    if (plot_file != "")
    {
        Eigen::VectorXd interpolantT = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
        Eigen::VectorXd interpolantE = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
        Eigen::VectorXd exact = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
        Eigen::VectorXd potential = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());
        Eigen::VectorXd elliptic = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());

        for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++)
        {
            exact(iV) = sol.value(mesh_ptr->vertex(iV)->coords());

            auto xV = mesh_ptr->vertex(iV)->coords();
            auto cList = mesh_ptr->vertex(iV)->get_cells();
            for (size_t ilT = 0; ilT < cList.size(); ilT++)
            {
                size_t iT = cList[ilT]->global_index();

                Eigen::VectorXd local_potential = PT[iT] * m_hho.restr(approx_Uvec, iT);
                Eigen::VectorXd local_elliptic = PT[iT] * m_hho.restr(interp_Uvec, iT);

                for (size_t i = 0; i < m_hho.local_highorder_dofs(iT); i++)
                {
                    potential(iV) += local_potential(i) * m_hho.highorder_basis(iT)->value(i, xV);
                    elliptic(iV) += local_elliptic(i) * m_hho.highorder_basis(iT)->value(i, xV);
                }
                for (size_t i = 0; i < m_hho.local_cell_dofs(iT); i++)
                {
                    interpolantT(iV) += m_hho.restr(interp_Uvec, iT)(i) * m_hho.highorder_basis(iT)->value(i, xV);
                }
            }
            potential(iV) = potential(iV) / (cList.size());
            elliptic(iV) = elliptic(iV) / (cList.size());
            interpolantT(iV) = interpolantT(iV) / (cList.size());

            auto eList = mesh_ptr->vertex(iV)->get_edges();
            for (size_t ilE = 0; ilE < cList.size(); ilE++)
            {
                size_t iE = eList[ilE]->global_index();

                Eigen::VectorXd interp_E = interp_Uvec.segment(m_hho.total_cell_dofs() + m_hho.global_offset_E(iE), m_hho.local_edge_dofs(iE));

                for (size_t i = 0; i < m_hho.local_edge_dofs(iE); i++)
                {
                    double tval = ((mesh_ptr->vertex(iV)->coords() - eList[ilE]->parameterisation().value(eList[ilE]->parameterisation().tmin)).norm() < 1E-12 ? eList[ilE]->parameterisation().tmin : eList[ilE]->parameterisation().tmax);
                    interpolantE(iV) += interp_E(i) * m_hho.edge_basis(iE)->value(i, tval);
                }
            }
            interpolantE(iV) = interpolantE(iV) / (eList.size());
        }

        VtuWriter plotdata(mesh_ptr);

        plotdata.write_to_vtu("exact-" + plot_file + ".vtu", exact);
        plotdata.write_to_vtu("potential-" + plot_file + ".vtu", potential);
        plotdata.write_to_vtu("elliptic-" + plot_file + ".vtu", elliptic);
        plotdata.write_to_vtu("IT-" + plot_file + ".vtu", interpolantT);
        plotdata.write_to_vtu("IE-" + plot_file + ".vtu", interpolantE);
        plotdata.write_to_vtu("diff-" + plot_file + ".vtu", potential - exact);
    }
}

void Model::local_poisson_operator(const size_t iT)
{
    Cell *cell = mesh_ptr->cell(iT);
    const size_t n_local_edges = cell->n_edges();
    const size_t n_local_cell_dofs = m_hho.local_cell_dofs(iT);
    const size_t n_local_highorder_dofs = m_hho.local_highorder_dofs(iT);
    const size_t n_local_bdry_dofs = m_hho.local_boundary_dofs(iT);
    const size_t n_local_dofs = n_local_cell_dofs + n_local_bdry_dofs;

    Eigen::MatrixXd MTT = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_highorder_dofs);

    QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

    for (size_t i = 0; i < n_local_cell_dofs; ++i)
    {
        for (size_t j = 0; j < n_local_highorder_dofs; ++j)
        {
            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                MTT(i, j) += quadT[iqn].w * m_hho.highorder_basis(iT)->value(i, quadT[iqn].x) * m_hho.highorder_basis(iT)->value(j, quadT[iqn].x);
            }
            // if ((i < j) && (j < n_local_cell_dofs))
            // {
            //     MTT(j, i) = MTT(i, j);
            // }
        }
    }

    // Calculate stiffness matrix
    Eigen::MatrixXd ST = Eigen::MatrixXd::Zero(n_local_highorder_dofs, n_local_highorder_dofs);

    for (size_t i = 0; i < n_local_highorder_dofs; ++i)
    {
        for (size_t j = 0; j < n_local_highorder_dofs; ++j)
        {
            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                ST(i, j) += quadT[iqn].w * m_hho.highorder_basis(iT)->derivative(i, quadT[iqn].x) * m_diffusion(cell) * (m_hho.highorder_basis(iT)->derivative(j, quadT[iqn].x)).transpose();
            }
            // ST(j, i) = ST(i, j);
        }
    }

    // Get L_T matrix to impose average condition on P_T later
    Eigen::VectorXd LT = (MTT.row(0)).transpose();
    Eigen::MatrixXd LTtLT = LT * (LT.transpose());
    double scalT = ST.norm() / LTtLT.norm();

    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_BDRY;
    std::vector<Eigen::Triplet<double>> triplets_K_M_BDRY_BDRY;
    // std::vector<Eigen::Triplet<double>> triplets_BDRY_SCALE_MAT;
    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_T;

    Eigen::SparseMatrix<double> M_BDRY_BDRY(n_local_bdry_dofs, n_local_bdry_dofs);
    Eigen::SparseMatrix<double> K_M_BDRY_BDRY(n_local_bdry_dofs, n_local_bdry_dofs);
    // Eigen::SparseMatrix<double> BDRY_SCALE_MAT(n_local_bdry_dofs, n_local_bdry_dofs);
    Eigen::SparseMatrix<double> M_BDRY_T(n_local_bdry_dofs, n_local_highorder_dofs);

    // Set cell - cell term of RHS of P_T
    Eigen::MatrixXd BP = Eigen::MatrixXd::Zero(n_local_highorder_dofs, n_local_dofs);
    BP.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs) = ST.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs) + scalT * LTtLT.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs); // cell bases contained in extended bases

    // double k_max = std::max(m_diffusion(cell)(0, 0), m_diffusion(cell)(1, 1));

    for (size_t iTF = 0; iTF < n_local_edges; iTF++)
    {
        Curve edge_param = cell->edge(iTF)->parameterisation();
        size_t iF = cell->edge(iTF)->global_index();

        const size_t n_local_edge_dofs = m_hho.local_edge_dofs(iF);
        const size_t offset_F = m_hho.local_offset_E(iT, iTF);
        const size_t offset_TF = n_local_cell_dofs + offset_F;

        QuadratureRule<double> quadF = m_hho.quadE(iF);

        for (size_t i = 0; i < n_local_edge_dofs; ++i)
        {
            for (size_t j = 0; j < n_local_edge_dofs; ++j)
            {
                double MFF_i_j = 0.0;
                double K_MFF_i_j = 0.0;
                for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
                {
                    Eigen::Vector2d normal = cell->edge_normal(iTF, quadF[iqn].x);
                    double temp = quadF[iqn].w * m_hho.edge_basis(iF)->value(i, quadF[iqn].x) * m_hho.edge_basis(iF)->value(j, quadF[iqn].x);
                    MFF_i_j += temp;
                    K_MFF_i_j += (m_diffusion(cell) * normal).dot(normal) * temp;
                    // K_MFF_i_j += k_max * temp;
                }
                triplets_M_BDRY_BDRY.emplace_back(offset_F + i, offset_F + j, MFF_i_j);
                triplets_K_M_BDRY_BDRY.emplace_back(offset_F + i, offset_F + j, K_MFF_i_j);
                // if (i < j)
                // {
                //     triplets_M_BDRY_BDRY.emplace_back(offset_F + j, offset_F + i, MFF_i_j);
                // }
            }
            for (size_t l = 0; l < n_local_highorder_dofs; ++l)
            {
                double MFT_i_l = 0.0;
                for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
                {
                    MFT_i_l += quadF[iqn].w * m_hho.edge_basis(iF)->value(i, quadF[iqn].x) * m_hho.highorder_basis(iT)->value(l, edge_param.value(quadF[iqn].x));
                }
                triplets_M_BDRY_T.emplace_back(offset_F + i, l, MFT_i_l);
            }
        }

        for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
        {
            for (size_t i = 0; i < n_local_highorder_dofs; ++i)
            {
                double grad_w_dot_nTF = m_hho.highorder_basis(iT)->derivative(i, edge_param.value(quadF[iqn].x)) * m_diffusion(cell) * (cell->edge_normal(iTF, quadF[iqn].x));
                for (size_t j = 0; j < n_local_cell_dofs; ++j)
                {
                    BP(i, j) -= quadF[iqn].w * grad_w_dot_nTF * m_hho.highorder_basis(iT)->value(j, edge_param.value(quadF[iqn].x));
                }
                for (size_t l = 0; l < n_local_edge_dofs; ++l)
                {
                    BP(i, offset_TF + l) += quadF[iqn].w * grad_w_dot_nTF * m_hho.edge_basis(iF)->value(l, quadF[iqn].x);
                }
            }
        }
    }

    M_BDRY_BDRY.setFromTriplets(std::begin(triplets_M_BDRY_BDRY), std::end(triplets_M_BDRY_BDRY));
    K_M_BDRY_BDRY.setFromTriplets(std::begin(triplets_K_M_BDRY_BDRY), std::end(triplets_K_M_BDRY_BDRY));
    M_BDRY_T.setFromTriplets(std::begin(triplets_M_BDRY_T), std::end(triplets_M_BDRY_T));

    PT[iT] = (ST + scalT * LTtLT).ldlt().solve(BP);
    Eigen::MatrixXd ATCONS = PT[iT].transpose() * ST * PT[iT];

    Eigen::MatrixXd DT = MTT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs).inverse() * MTT * PT[iT];

    // Eigen::MatrixXd DT = MTT * PT[iT];
    DT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) -= Eigen::MatrixXd::Identity(n_local_cell_dofs, n_local_cell_dofs);

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    // Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_BDRY_BDRY);

    // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    // solver.analyzePattern(M_BDRY_BDRY);
    // solver.factorize(M_BDRY_BDRY);

    Eigen::MatrixXd M_BDRY_BDRY_INV_M_BDRY_T = solver.solve(M_BDRY_T);
    // Eigen::MatrixXd M_BDRY_BDRY_INV_M_BDRY_T = M_BDRY_T;
    Eigen::MatrixXd STAB_OPERATOR = M_BDRY_BDRY_INV_M_BDRY_T.topLeftCorner(n_local_bdry_dofs, n_local_highorder_dofs) * PT[iT];
    STAB_OPERATOR.bottomRightCorner(n_local_bdry_dofs, n_local_bdry_dofs) -= Eigen::MatrixXd::Identity(n_local_bdry_dofs, n_local_bdry_dofs);

    Eigen::MatrixXd ATSTAB = DT.transpose() * ST.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) * DT + (1.0 / cell->diam()) * STAB_OPERATOR.transpose() * K_M_BDRY_BDRY * STAB_OPERATOR;
    // Eigen::MatrixXd ATSTAB = (1.0 / cell->diam()) * STAB_OPERATOR.transpose() * K_M_BDRY_BDRY * STAB_OPERATOR;

    // Eigen::MatrixXd ATSTAB = DT.transpose() * ST.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) * DT + (1.0 / cell->diam()) * STAB_OPERATOR.transpose() * BDRY_SCALE_MAT * STAB_OPERATOR;

    AT[iT] = ATCONS + ATSTAB;
}

void Model::local_source_term(const size_t iT, Eigen::VectorXd &bT)
{
    bT = Eigen::VectorXd::Zero(m_hho.local_cell_dofs(iT) + m_hho.local_boundary_dofs(iT));

    QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);
    for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
    {
        // double source_weight = quadT[iqn].w * m_src.value(quadT[iqn].x);
        double source_weight = quadT[iqn].w * m_src(quadT[iqn].x, mesh_ptr->cell(iT));
        for (size_t i = 0; i < m_hho.local_cell_dofs(iT); ++i)
        {
            bT(i) += source_weight * m_hho.highorder_basis(iT)->value(i, quadT[iqn].x);
        }
    }
}