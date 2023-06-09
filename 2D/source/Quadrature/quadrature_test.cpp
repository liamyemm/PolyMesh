// internal libraries
#include "QuadratureRule.hpp"
#include "basis.hpp"
#include "function.hpp"
#include "MeshBuilder2D.hpp"
#include "MeshReaderTyp2.hpp"
#include "mesh_convert.hpp"

#include "QuadHandler.hpp"
#include "generate_quadrature_rule.hpp"

// boost libraries
#include <boost/program_options.hpp> // boost::program_options
#include <boost/timer/timer.hpp>     // boost::cpu_timer

// std libraries
#include <iostream> // std::cout
#include <iomanip>  // std::setprecision

using namespace PolyMesh2D::Quadrature;

int main(const int argc, const char **argv)
{
    // QuadratureRule<double> quad(gauss_jacobi(0, -0.95, 300));

    // std::function<double(double)> test_g = [](double x) -> double
    // {
    //     return std::pow(x, 3) - 17.0 * std::pow(x, 8) + 14.0 * x - 12.0;
    // };

    // std::function<double(double)> test_f = [](double x) -> double
    // {
    //     return std::pow(x, -0.95) * (std::pow(x, 3) - 17.0 * std::pow(x, 8) + 14.0 * x - 12.0);
    // };

    // double integral1 = 0.0;

    // for(size_t qn = 0; qn < quad.size(); ++qn)
    // {
    //     integral1 += quad.weight(qn) * test_g(quad.point(qn));
    // }

    // normalise_weights(quad, 0, -0.95);

    // double integral2 = 0.0;

    // for(size_t qn = 0; qn < quad.size(); ++qn)
    // {
    //     integral2 += quad.weight(qn) * test_f(quad.point(qn));
    // }

    // std::cout << "\n" << std::abs(integral1 - -228.4505990564437) << "\n";

    // std::cout << "\n" << std::abs(integral2 - -228.4505990564437) << "\n";

    // Build mesh and reorder edges
    PolyMesh2D::StraightMesh::MeshBuilder builder = PolyMesh2D::StraightMesh::MeshBuilder("/home/liam/github/codes/liamyemm/PolyMesh/2D/typ2_meshes/mesh2_1.typ2");
    std::unique_ptr<PolyMesh2D::StraightMesh::Mesh> straight_mesh;

    try
    {
        straight_mesh = builder.build_the_mesh();
    }
    catch (std::string msg)
    {
        std::cerr << msg;
    }

    std::unique_ptr<PolyMesh2D::CurvedMesh::Mesh> curved_mesh = PolyMesh2D::MeshTransform::mesh_convert(straight_mesh.get());

    straight_mesh.reset(); // delete straight mesh

    PolyMesh2D::Quadrature::QuadHandler<PolyMesh2D::CurvedMesh::Mesh> quad_handle(curved_mesh.get(), 18, 18);

    // double scalar_pole_order1 = 1.5;
    double scalar_pole_order1 = -0.5;
    double scalar_pole_order2 = 0.5;

    // double vector_pole_order1 = -0.5;
    double vector_pole_order1 = 1.4;
    double vector_pole_order2 = 1.9;

    std::function<double(Eigen::Vector2d)> scalar_pole_1_lambda = [&scalar_pole_order1](Eigen::Vector2d x) -> double
    {
        return std::pow(x.norm(), -scalar_pole_order1);
    };

    std::function<Eigen::RowVector2d(Eigen::Vector2d)> D_scalar_pole_1_lambda = [&scalar_pole_order1](Eigen::Vector2d x) -> Eigen::RowVector2d
    {
        return 0.5 * std::pow(x.norm(), -scalar_pole_order1 - 1.0) * x.transpose() / x.norm();
    };

    std::function<double(Eigen::Vector2d)> scalar_2_pole_lambda = [&scalar_pole_order2](Eigen::Vector2d x) -> double
    {
        return std::pow((x - Eigen::Vector2d(1, 1)).norm(), -scalar_pole_order2);
    };

    std::function<Eigen::Vector2d(Eigen::Vector2d)> vector_pole_1_lambda = [&vector_pole_order1](Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return Eigen::Vector2d(1, 1) * std::pow(x.norm(), -vector_pole_order1);
    };

    std::function<Eigen::Vector2d(Eigen::Vector2d)> vector_pole_2_lambda = [&vector_pole_order2](Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return Eigen::Vector2d(1, 1) * std::pow((x - Eigen::Vector2d(1, 1)).norm(), -vector_pole_order2);
    };

    std::function<double(Eigen::Vector2d)> polynomial_lambda = [](Eigen::Vector2d x) -> double
    {
        return (std::pow(x(0) * x(1), 3) - 17.0 * std::pow(x(0), 8) + 12.0 * std::pow(x(1), 6) * x(0) + 14.0 * x(0) * x(1) - 12.0);
    };

    std::function<Eigen::RowVector2d(Eigen::Vector2d)> D_polynomial_lambda = [](Eigen::Vector2d x) -> Eigen::RowVector2d
    {
        return Eigen::RowVector2d(-136.0 * std::pow(x(0), 7) + 3.0 * std::pow(x(0) * x(1), 2) * x(1) + 2 * x(1) * (7.0 + 6.0 * std::pow(x(1), 5)), x(0) * (14.0 + 3.0 * std::pow(x(0) * x(1), 2) + 72.0 * std::pow(x(1), 5)));
    };

    PolyMesh2D::Functional::Function<2, 1> scalar_pole1(scalar_pole_1_lambda, D_scalar_pole_1_lambda);
    PolyMesh2D::Functional::Function<2, 1> scalar_pole2(scalar_2_pole_lambda);
    PolyMesh2D::Functional::Function<2, 2> vector_pole1(vector_pole_1_lambda);
    PolyMesh2D::Functional::Function<2, 2> vector_pole2(vector_pole_2_lambda);
    PolyMesh2D::Functional::Function<2, 1> polynomial(polynomial_lambda, D_polynomial_lambda);

    PolyMesh2D::Functional::Pole<Eigen::Vector2d> scalar_pole_1;
    scalar_pole_1.location = Eigen::Vector2d::Zero();
    scalar_pole_1.order = scalar_pole_order1;
    scalar_pole1.add_pole(scalar_pole_1);

    PolyMesh2D::Functional::Pole<Eigen::Vector2d> scalar_pole_2;
    scalar_pole_2.location = Eigen::Vector2d(1, 1);
    scalar_pole_2.order = scalar_pole_order2;
    scalar_pole2.add_pole(scalar_pole_2);

    PolyMesh2D::Functional::Pole<Eigen::Vector2d> vector_pole_1;
    vector_pole_1.location = Eigen::Vector2d::Zero();
    vector_pole_1.order = vector_pole_order1;
    vector_pole1.add_pole(vector_pole_1);

    PolyMesh2D::Functional::Pole<Eigen::Vector2d> vector_pole_2;
    vector_pole_2.location = Eigen::Vector2d(1, 1);
    vector_pole_2.order = vector_pole_order2;
    vector_pole2.add_pole(vector_pole_2);

    double integralScalarPole1 = 0.0;
    double integralScalarPole2 = 0.0;
    double integralScalarPole1TimesPole2 = 0.0;
    double integralScalarPole1TimesPoly = 0.0;
    double integralScalarPole2TimesPoly = 0.0;
    double integralScalarPole1TimesPole2TimesPoly = 0.0;
    double integralDotVectorPole1VectorPole2 = 0.0;
    double integralDotVectorPole1VectorPole2TimesPoly = 0.0;

    // double int2minusint1 = 0.0;
    // double int2plusint1 = 0.0;

    for (auto &cell : curved_mesh->get_cells())
    {
        integralScalarPole1 += quad_handle.integrate(scalar_pole1, cell);
        integralScalarPole2 += quad_handle.integrate(scalar_pole2, cell);
        integralScalarPole1TimesPole2 += quad_handle.integrate(scalar_pole1 * scalar_pole2, cell);
        integralScalarPole1TimesPoly += quad_handle.integrate(scalar_pole1 * polynomial, cell);
        integralScalarPole2TimesPoly += quad_handle.integrate(scalar_pole2 * polynomial, cell);
        integralScalarPole1TimesPole2TimesPoly += quad_handle.integrate(scalar_pole1 * scalar_pole2 * polynomial, cell);
        integralDotVectorPole1VectorPole2 += quad_handle.integrate(PolyMesh2D::Functional::scalar_product(vector_pole1, vector_pole2), cell);
        integralDotVectorPole1VectorPole2TimesPoly += quad_handle.integrate(PolyMesh2D::Functional::scalar_product(vector_pole1, vector_pole2) * polynomial, cell);

        // int2minusint1 += quad_handle.integrate((PolyMesh2D::Functional::scalar_product(vector_pole1, vector_pole2) * polynomial) - (PolyMesh2D::Functional::scalar_product(vector_pole1, vector_pole2)), cell);
        // int2plusint1 += quad_handle.integrate((PolyMesh2D::Functional::scalar_product(vector_pole1, vector_pole2) * polynomial) + (PolyMesh2D::Functional::scalar_product(vector_pole1, vector_pole2)), cell);
    }

    std::cout << "Testing cell integration\n" << std::endl;

    // Output the results
    std::cout << "Integration error of scalarPole1: " << std::abs(integralScalarPole1 - 3.323584864723750) << std::endl;
    std::cout << "Integration error of scalarPole2: " << std::abs(integralScalarPole2 - 1.249986334329248) << std::endl;
    std::cout << "Integration error of scalarPole1 times scalarPole2: " << std::abs(integralScalarPole1TimesPole2 - 3.333970240854988) << std::endl; 
    std::cout << "Integration error of scalarPole1 times poly: " << std::abs(integralScalarPole1TimesPoly - -36.64098524761831) << std::endl;
    std::cout << "Integration error of scalarPole2 times poly: " << std::abs(integralScalarPole2TimesPoly - -10.783391135827567) << std::endl;
    std::cout << "Integration error of scalarPole1 times scalarPole2 times poly: " << std::abs(integralScalarPole1TimesPole2TimesPoly - -35.34406818583255) << std::endl;
    std::cout << "Integration error of vectorPole1 dot vectorPole2 (int1): " << std::abs(integralDotVectorPole1VectorPole2 - 36.45859882562320) << std::endl;
    std::cout << "Integration error of vectorPole1 dot vectorPole2 times poly (int2): " << std::abs(integralDotVectorPole1VectorPole2TimesPoly - -94.92387830666720) << std::endl;

    // std::cout << "Integration error of int2minusint1: " << std::abs(int2minusint1 - (-94.92387830666720 - 36.45859882562320)) << std::endl;
    // std::cout << "Integration error of int2plusint1: " << std::abs(int2plusint1 - (-94.92387830666720 + 36.45859882562320)) << std::endl;

    std::cout << "Testing edge integration\n" << std::endl;

    std::cout << curved_mesh->edge(0)->vertex(0)->coords() << "\n\n" << curved_mesh->edge(0)->vertex(1)->coords() << "\n";
    std::cout << curved_mesh->edge(1)->vertex(0)->coords() << "\n\n" << curved_mesh->edge(1)->vertex(1)->coords() << "\n";

    auto trace = PolyMesh2D::Functional::trace(scalar_pole1, curved_mesh->edge(0)->parameterisation());
    std::cout << "\n" << std::abs(quad_handle.integrate(trace, curved_mesh->edge(0)) - 1.0 / 12.0) << "\n";

    auto neumann_trace = PolyMesh2D::Functional::neumann_trace(scalar_pole1 * polynomial, curved_mesh->edge(0)->parameterisation());
    std::cout << "\n" << std::abs(quad_handle.integrate(neumann_trace, curved_mesh->edge(0)) - 717.0 / 4096.0) << "\n";

    std::cout << "\n" << curved_mesh->edge(0)->parameterisation().value(curved_mesh->edge(0)->parameterisation().tmin) << "\n";
    std::cout << "\n" << curved_mesh->edge(1)->parameterisation().value(curved_mesh->edge(0)->parameterisation().tmin) << "\n";

    auto trace2 = PolyMesh2D::Functional::trace(scalar_pole1, curved_mesh->edge(1)->parameterisation());
    std::cout << "\n" << std::abs(quad_handle.integrate(trace2, curved_mesh->edge(1)) - 1.0 / 12.0) << "\n";

    auto neumann_trace2 = PolyMesh2D::Functional::neumann_trace(scalar_pole1 * polynomial, curved_mesh->edge(1)->parameterisation());
    std::cout << "\n" << std::abs(quad_handle.integrate(neumann_trace2, curved_mesh->edge(1)) - 7.0 / 40.0) << "\n";

    auto new_func = PolyMesh2D::Functional::scalar_product(vector_pole1, PolyMesh2D::Functional::gradient(scalar_pole1 * polynomial));
    // auto new_func = PolyMesh2D::Functional::gradient(scalar_pole1);
    // auto new_func = PolyMesh2D::Functional::gradient(polynomial);

    double new_int = 0.0;
    // Eigen::Vector2d vec_int = Eigen::Vector2d::Zero();

    // std::function<double(Eigen::Vector2d)> grad = [scalar_pole1](const Eigen::Vector2d &x) -> double
    // {
    //     return scalar_pole1.derivative(x)(0);
    // };

    // PolyMesh2D::Functional::Function<2, 1> d0(grad);

    // PolyMesh2D::Functional::Pole<Eigen::Vector2d> grad_pole(Eigen::Vector2d::Zero(), scalar_pole_order1 + 1.0);

    // d0.add_pole(grad_pole);

    for (auto &cell : curved_mesh->get_cells())
    {
        new_int += quad_handle.integrate(new_func, cell);
        // vec_int += quad_handle.integrate(new_func, cell);
        // new_int += quad_handle.integrate(d0, cell);
    }
    std::cout << "Integration error of new_int: " << std::abs(new_int - -110.2030975449566) << std::endl;
    // std::cout << "Integration error of new_int: " << std::abs(new_int - -7.80851) << std::endl;
    // std::cout << "Integration error of new_int: " << (vec_int - Eigen::Vector2d(0.4031034820621265, 0.4031034820621265)).norm() << std::endl;
    // std::cout << "Integration error of new_int: " << std::abs(new_int - 0.4031034820621265) << std::endl;
    // std::cout << "Integration error of new_int: " << (vec_int - Eigen::Vector2d(-8.035714285714286, 13.25000000000000)).norm() << std::endl;
    // std::cout << "Integration error of new_int: " << std::abs(new_int - 13.25000000000000) << std::endl;

    // std::cout << "\n" << std::abs(mesh_integral - -9.469246031746032) << "\n";
    // std::cout << "\n" << std::abs(mesh_integral - 405.0640705787629) << "\n";
    // std::cout << "\n" << std::abs(mesh_integral - 3.323584864723750) << "\n";
    // std::cout << "\n" << std::abs(mesh_integral - -36.64098524761831) << "\n";

    std::cout << curved_mesh->edge(13)->vertex(0)->coords() << "\n\n\n\n" << curved_mesh->edge(13)->vertex(1)->coords() << "\n\n\n";

    auto vF_dot_n(PolyMesh2D::Functional::dot_n(PolyMesh2D::Functional::trace(vector_pole1, curved_mesh->edge(0)->parameterisation()), curved_mesh->edge(0)->parameterisation()));
    auto pT_vF_dot_n(vF_dot_n * PolyMesh2D::Functional::trace(scalar_pole1 * polynomial, curved_mesh->edge(0)->parameterisation()));
    std::cout << "\n" << std::abs(quad_handle.integrate(pT_vF_dot_n, curved_mesh->edge(0)) - -104.4660675955349) << "\n";

    auto vF_dot_n2(PolyMesh2D::Functional::dot_n(PolyMesh2D::Functional::trace(vector_pole1, curved_mesh->edge(13)->parameterisation()), curved_mesh->edge(13)->parameterisation()));
    auto pT_vF_dot_n2(vF_dot_n2 * PolyMesh2D::Functional::trace(scalar_pole1 * polynomial, curved_mesh->edge(13)->parameterisation()));
    std::cout << "\n" << std::abs(quad_handle.integrate(pT_vF_dot_n2, curved_mesh->edge(13)) - -7.497891388881993) << "\n";
}