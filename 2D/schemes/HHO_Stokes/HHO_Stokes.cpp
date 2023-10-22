#include "HHO_Stokes.hpp"

// internal libraries
#include "parallel_for.hpp"
#include "vtu_writer.hpp"
#include "QuadratureRule.hpp"
#include "basis.hpp"
#include "function.hpp"
#include "MeshBuilder2D.hpp"
#include "MeshReaderTyp2.hpp"

#include "mesh_convert.hpp"
#include "mesh_cutter.hpp"
#include "../HHO_Poisson/TestCase.hpp"

// boost libraries
#include <boost/program_options.hpp> // boost::program_options
#include <boost/timer/timer.hpp>     // boost::cpu_timer

// std libraries
#include <iostream> // std::cout
#include <iomanip>  // std::setprecision

#include <numeric> // std::accumulate

#include "vtk_writer.hpp"

using namespace PolyMesh2D::HHOSTOKES;

Eigen::Matrix2d tensor_product(const Eigen::Vector2d &v1, const Eigen::Vector2d &v2)
{
    Eigen::Matrix2d tensorProduct;
    tensorProduct << v1(0) * v2(0), v1(0) * v2(1), v1(1) * v2(0), v1(1) * v2(1);
    return tensorProduct;
}

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

    std::unique_ptr<PolyMesh2D::CurvedMesh::Mesh> curved_mesh = PolyMesh2D::MeshTransform::mesh_convert(straight_mesh.get());

    straight_mesh.reset(); // delete straight mesh

    std::vector<PolyMesh2D::Circle> circles;

    {
        // PolyMesh2D::Circle circle1(Eigen::Vector2d(0.4, 0.275), 0.01);
        // PolyMesh2D::Circle circle2(Eigen::Vector2d(0.2, 0.25), 0.02);
        // PolyMesh2D::Circle circle3(Eigen::Vector2d(0.7, 0.5), 0.075);
        // PolyMesh2D::Circle circle4(Eigen::Vector2d(0.3, 0.75), 0.015);

        // circles = {circle1, circle2, circle3, circle4};

        PolyMesh2D::Circle circle(Eigen::Vector2d(0.5, 0.5), params.radius);
        circles = {circle};
    }

    for (auto &circle : circles)
    {
        PolyMesh2D::MeshCutter mesh_cutter(curved_mesh.get(), circle.level_set(), circle.param(), false, true);
        mesh_cutter.cut_mesh();
    }

    std::ofstream mesh_out("mesh_plot.dat");
    curved_mesh->plot_mesh(&mesh_out, 20);
    mesh_out.close();

    assert(curved_mesh->test());
    exit(1);

    // Reordering
    PolyMesh2D::HHOPOISSON::reorder_edges(curved_mesh.get());

    // Print the data
    std::cout << "[Scheme] Data:\n";
    std::cout << "     No. cells = " << curved_mesh->n_cells() << ", No. edges = " << curved_mesh->n_edges() << ", No. vertices = " << curved_mesh->n_vertices() << "\n";
    // std::cout << "     TestCase = " << params.test_case << "\n";
    std::cout << "     Mesh = " << params.mesh_name << "\n";
    std::cout << "     Degrees: edge = " << params.edge_degree << "; cell = " << params.cell_degree << "\n";
    std::cout << "     Using threads = " << (params.use_threads ? "true" : "false") << std::endl;

    PolyMesh2D::StokesCore stokes(curved_mesh.get(), params.cell_degree, params.edge_degree, params.use_threads, params.orthonormalise);

    // PolyMesh2D::StokesTestSquare stokes_test_square;

    // VectorFunction2D source(stokes_test_square.source());
    // VectorFunction2D u_func(stokes_test.u());
    // VectorFunction2D laplace_u_func(stokes_test.laplace_u());
    // ScalarFunction2D p_func(stokes_test.p());

    std::function<Eigen::Vector2d(Eigen::Vector2d)> src = [](const Eigen::Vector2d x) -> Eigen::Vector2d
    {
        return Eigen::Vector2d::Zero();
    };

    VectorFunction2D source(src);

    // source.add_pole(PolyMesh2D::Functional::Pole<Eigen::Vector2d>(Eigen::Vector2d::Zero(), 2.0 - lambda));

    // VectorFunction2D u_func(zero_u, D_zero_u);
    // VectorFunction2D laplace_u_func(zero_u, D_zero_u);

    // VectorFunction2D u_func(circle.u() + stokes_test_square.u());
    // VectorFunction2D laplace_u_func(circle.laplace_u() + stokes_test_square.laplace_u());
    // ScalarFunction2D p_func(circle.p(0.0) + stokes_test_square.p(0.0));

    VectorFunction2D u_func(circles[0].u());
    ScalarFunction2D p_func(circles[0].p(0.0));

    // VectorFunction2D u_func(circle.u());
    // VectorFunction2D laplace_u_func(circle.laplace_u());
    // ScalarFunction2D p_func(circle.p(0.0));

    // std::function<double(Eigen::Vector2d)> test = [sing, stokes_test](const Eigen::Vector2d x) -> double
    // {
    //     return (-sing.laplace_u().value(x) + sing.p(0.0).derivative(x).transpose()).squaredNorm();
    // };
    // ScalarFunction2D test_func(test);

    // ScalarFunction2D test_func(PolyMesh2D::Functional::curl(sing.laplace_u()));

    // double integral = 0.0;
    // auto handle = stokes.get_quad_handle();
    // for(auto & cell : curved_mesh->get_cells())
    // {
    //     integral += handle.integrate(p_func, cell);
    // }

    // std::cout << "\n" << std::setprecision(16) << integral << "\n";
    // exit(1);

    double cutoff = 0.1;

#if 1
    for (auto &circle : circles)
    {
        for (size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            if (circle.enrich_cell(curved_mesh->cell(iT), cutoff))
            {
                stokes.enrich_highorder_basis(iT, circle.u());
                stokes.enrich_pressure_basis(iT, circle.p(0.0));
                stokes.enrich_cell_basis(iT, circle.laplace_u());
            }
        }

        for (size_t iE = 0; iE < curved_mesh->n_edges(); ++iE)
        {
            if (circle.enrich_edge(curved_mesh->edge(iE), cutoff))
            {
                // if((circle.u_linearly_independent(curved_mesh->edge(iE))) || (params.edge_degree == 0))
                // {
                stokes.enrich_edge_basis(iE, PolyMesh2D::Functional::neumann_trace(circle.u(), curved_mesh->edge(iE)->parameterisation()));
                // }

                // if(circle.p_linearly_independent(curved_mesh->edge(iE)))
                // {
                stokes.enrich_edge_basis(iE, PolyMesh2D::Functional::times_n(PolyMesh2D::Functional::trace(circle.p(0.0), curved_mesh->edge(iE)->parameterisation()), curved_mesh->edge(iE)->parameterisation()));
                // }
            }
        }
    }
#endif

    Model model(stokes, source);

    model.assemble(params.use_threads);

    Eigen::VectorXd interp(stokes.interpolate(u_func, p_func));

    size_t bdry_vector_dofs = 0;
    for (size_t ibF = 0; ibF < curved_mesh->n_b_edges(); ++ibF)
    {
        const size_t iF = curved_mesh->b_edge(ibF)->global_index();
        bdry_vector_dofs += stokes.local_edge_dofs(iF);
    }

    // set boundary dofs
    // Eigen::VectorXd UDir = interp.segment(stokes.total_cell_dofs() + stokes.total_edge_dofs() - bdry_vector_dofs, bdry_vector_dofs);

    Eigen::VectorXd UDir = Eigen::VectorXd::Zero(bdry_vector_dofs);

    std::function<Eigen::Vector2d(double)> dirichlet_bc = [](const double x) -> Eigen::Vector2d
    {
        return Eigen::Vector2d(1.0, 0.0);
    };

    PolyMesh2D::Functional::VectorFunction1D BC(dirichlet_bc);

    for (auto &b_edge : curved_mesh->get_b_edges())
    {
        if (b_edge->is_straight())
        {
            UDir.segment(stokes.global_offset_E(b_edge->global_index()) - stokes.total_edge_dofs() + bdry_vector_dofs, stokes.local_edge_dofs(b_edge->global_index())) = stokes.l2_edge_projection(b_edge->global_index(), BC);
        }
        else
        {
            UDir.segment(stokes.global_offset_E(b_edge->global_index()) - stokes.total_edge_dofs() + bdry_vector_dofs, stokes.local_edge_dofs(b_edge->global_index())) = Eigen::VectorXd::Zero(stokes.local_edge_dofs(b_edge->global_index()));
        }
    }

    Eigen::VectorXd sol(model.solve(UDir));

    std::vector<Eigen::VectorXd> elliptic_projectors;
    std::vector<Eigen::VectorXd> l2_velocity_projectors;
    elliptic_projectors.reserve(curved_mesh->n_cells());
    l2_velocity_projectors.reserve(curved_mesh->n_cells());
    for (size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
    {
        elliptic_projectors.push_back(stokes.elliptic_projection(iT, u_func));
        l2_velocity_projectors.push_back(stokes.l2_highorder_projection(iT, u_func));
    }

    auto errors = model.compute_errors(sol, interp, elliptic_projectors, l2_velocity_projectors);

    double l2_error = errors[0];
    double h1_error = errors[1];
    double energy_error = errors[2];
    double pressure_error = errors[3];
    double divergence_error = errors[4];

    std::cout << "     L2 Error = " << l2_error << "\n";
    std::cout << "     H1 Error = " << h1_error << "\n";
    std::cout << "     Energy Error = " << energy_error << "\n";
    std::cout << "     Pressure Error = " << pressure_error << "\n";
    std::cout << "     Divergence Error = " << divergence_error << "\n";

    model.plot(sol);

    std::cout << "\n[Scheme] Writing solution to file\n";

    double global_dofs = stokes.total_edge_dofs() - bdry_vector_dofs + curved_mesh->n_cells();

    std::ofstream out("results.txt");
    out << "EdgeDegree: " << params.edge_degree << "\n";
    out << "CellDegree: " << params.cell_degree << "\n";
    out << "L2Error: " << l2_error << "\n";
    out << "H1Error: " << h1_error << "\n";
    out << "EnergyError: " << energy_error << "\n";
    out << "PressureError: " << pressure_error << "\n";
    out << "MeshSize: " << curved_mesh->h_max() << "\n";
    out << "NbCells: " << curved_mesh->n_cells() << "\n";
    out << "NbEdges: " << curved_mesh->n_edges() << "\n";
    out << "NbInternalEdges: " << curved_mesh->n_edges() - curved_mesh->n_b_edges() << "\n";
    out << "DOFs: " << global_dofs << "\n";
    out << "Radius: " << params.radius << "\n";
    out.close();

    return 0;
}

ModelParameters::ModelParameters(const int argc, const char **argv)
{
    namespace po = boost::program_options;

    // Program options
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Produce help message")("mesh,m", po::value<std::string>(), "Set the mesh")("celldegree,l", po::value<unsigned>(), "Set the degree of the cell polynomials")("edgedegree,k", po::value<unsigned>(), "Set the degree of the edge polynomials")("plot,p", po::value<std::string>(), "Plot to file")("use_threads,u", po::value<bool>(), "Using multithreading")("orthonormalise,o", po::value<bool>(), "Orthonormalise the basis functions")("radius,r", po::value<double>(), "Set the circle radius");

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
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "../../../typ2_meshes/hexa1_3.typ2");
    // mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "../../../typ2_meshes/mesh2_2_transformed.typ2");

    // Get polynomial degrees
    edge_degree = (vm.count("edgedegree") ? vm["edgedegree"].as<unsigned>() : 0);
    cell_degree = (vm.count("celldegree") ? vm["celldegree"].as<unsigned>() : edge_degree);

    // if ((std::abs((int)cell_degree - (int)edge_degree) > 1))
    // {
    //     std::cout << "Invalid choice of polynomial degrees\n";
    //     exit(1);
    // }

    // Get plot file
    plot_file = (vm.count("plot") ? vm["plot"].as<std::string>() : "");

    // Get use_threads
    use_threads = (vm.count("use_threads") ? vm["use_threads"].as<bool>() : true);

    orthonormalise = (vm.count("orthonormalise") ? vm["orthonormalise"].as<bool>() : true);

    radius = (vm.count("radius") ? vm["radius"].as<double>() : 0.1);
}

Model::Model(const StokesCore &stokes, const VectorFunction2D &src) : m_stokes(stokes), m_src(src), mesh_ptr(m_stokes.get_mesh()) {}

void Model::local_stokes_operator(const size_t iT, Eigen::MatrixXd &AT_loc, Eigen::MatrixXd &RT_loc)
{
    const size_t n_local_velocity_cell_dofs = m_stokes.local_cell_dofs(iT);
    const size_t n_local_velocity_bdry_dofs = m_stokes.local_boundary_dofs(iT);
    const size_t n_local_velocity_highorder_dofs = m_stokes.local_highorder_dofs(iT);
    const size_t n_local_pressure_dofs = m_stokes.local_pressure_dofs(iT);

    // QuadratureRule<Eigen::Vector2d> quadT(m_stokes.quadT(iT));

    CurvedMesh::Cell *cell = mesh_ptr->cell(iT);
    auto quad_handle(m_stokes.get_quad_handle());

    PolyMesh2D::Functional::ScalarFamily2D *pressure_basis(m_stokes.pressure_basis(iT));
    PolyMesh2D::Functional::VectorFamily2D *velocity_cell_basis(m_stokes.cell_basis(iT));
    PolyMesh2D::Functional::VectorFamily2D *velocity_highorder_basis(m_stokes.highorder_basis(iT));

    Eigen::MatrixXd DT_RHS = Eigen::MatrixXd::Zero(n_local_pressure_dofs, n_local_velocity_cell_dofs + n_local_velocity_bdry_dofs);
    Eigen::MatrixXd RT_RHS = Eigen::MatrixXd::Zero(n_local_velocity_highorder_dofs, n_local_velocity_cell_dofs + n_local_velocity_bdry_dofs);

    const bool sym = true;
    const bool not_sym = false;

    Eigen::MatrixXd M_CELL_CELL = quad_handle.l2_product(*velocity_cell_basis, *velocity_cell_basis, cell, sym);
    Eigen::MatrixXd M_CELL_HIGH = quad_handle.l2_product(*velocity_cell_basis, *velocity_highorder_basis, cell, not_sym);

    Eigen::MatrixXd ST = quad_handle.h1_product(*velocity_highorder_basis, *velocity_highorder_basis, cell, sym);

    for (size_t i = 0; i < n_local_velocity_cell_dofs; ++i)
    {
        for (size_t j = 0; j < n_local_pressure_dofs; ++j)
        {
            DT_RHS(j, i) = -quad_handle.integrate(Functional::scalar_product(velocity_cell_basis->ancestor().function(i), Functional::gradient(pressure_basis->ancestor().function(j))), cell);
        }
    }

    Eigen::MatrixXd tmp = DT_RHS.topLeftCorner(n_local_pressure_dofs, n_local_velocity_cell_dofs);
    DT_RHS.topLeftCorner(n_local_pressure_dofs, n_local_velocity_cell_dofs) = pressure_basis->matrix() * tmp * velocity_cell_basis->matrix().transpose();

    Eigen::MatrixXd LT = (M_CELL_HIGH.topRightCorner(2, n_local_velocity_highorder_dofs)).transpose(); // need both components of the constant term
    Eigen::MatrixXd LTtLT = LT * (LT.transpose());

    RT_RHS.topLeftCorner(n_local_velocity_highorder_dofs, n_local_velocity_cell_dofs) = quad_handle.h1_product(*velocity_highorder_basis, *velocity_cell_basis, cell, not_sym);

    Eigen::MatrixXd L_HIGH_CELL = (M_CELL_CELL.block(0, 0, n_local_velocity_cell_dofs, 2) * M_CELL_HIGH.block(0, 0, 2, n_local_velocity_highorder_dofs)).transpose();
    Eigen::MatrixXd L_HIGH_HIGH = M_CELL_HIGH.block(0, 0, 2, n_local_velocity_highorder_dofs).transpose() * M_CELL_HIGH.block(0, 0, 2, n_local_velocity_highorder_dofs);

    double scalT = ST.norm() / L_HIGH_HIGH.norm();

    RT_RHS.topLeftCorner(n_local_velocity_highorder_dofs, n_local_velocity_cell_dofs) += scalT * L_HIGH_CELL;

    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_BDRY;
    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_HIGH;

    Eigen::SparseMatrix<double> M_BDRY_BDRY(n_local_velocity_bdry_dofs, n_local_velocity_bdry_dofs);
    Eigen::SparseMatrix<double> M_BDRY_HIGH(n_local_velocity_bdry_dofs, n_local_velocity_highorder_dofs);

    for (size_t iTE = 0; iTE < mesh_ptr->cell(iT)->n_edges(); iTE++)
    {
        const size_t iE = mesh_ptr->cell(iT)->edge(iTE)->global_index();
        const size_t n_local_velocity_edge_dofs = m_stokes.local_edge_dofs(iE);
        PolyMesh2D::Functional::VectorFamily1D *velocity_edge_basis(m_stokes.edge_basis(iE));

        const size_t local_offset = n_local_velocity_cell_dofs + m_stokes.local_offset_E(iT, iTE);

        Curve edge_param = mesh_ptr->cell(iT)->edge(iTE)->parameterisation();

        { // scope to set DT_RHS
            Eigen::MatrixXd pT_vF_dot_n = Eigen::MatrixXd::Zero(n_local_pressure_dofs, n_local_velocity_edge_dofs);
            for (size_t j = 0; j < n_local_velocity_edge_dofs; ++j)
            {
                auto vF_dot_n(Functional::dot_n(velocity_edge_basis->ancestor().function(j), edge_param));
                for (size_t i = 0; i < n_local_pressure_dofs; ++i)
                {
                    pT_vF_dot_n(i, j) = quad_handle.integrate(vF_dot_n * Functional::trace(pressure_basis->ancestor().function(i), edge_param), cell->edge(iTE));
                }
            }
            DT_RHS.block(0, local_offset, n_local_pressure_dofs, n_local_velocity_edge_dofs) = cell->edge_orientation(iTE) * pressure_basis->matrix() * pT_vF_dot_n * velocity_edge_basis->matrix().transpose();
        }
        { // scope to set RT_RHS
            Eigen::MatrixXd gradw_nTF_vF = Eigen::MatrixXd::Zero(n_local_velocity_highorder_dofs, n_local_velocity_edge_dofs);
            Eigen::MatrixXd gradw_nTF_vT = Eigen::MatrixXd::Zero(n_local_velocity_highorder_dofs, n_local_velocity_cell_dofs);
            for (size_t i = 0; i < n_local_velocity_highorder_dofs; ++i)
            {
                auto gradw_nTF(Functional::neumann_trace(velocity_highorder_basis->ancestor().function(i), edge_param));
                for (size_t j = 0; j < n_local_velocity_edge_dofs; ++j)
                {
                    gradw_nTF_vF(i, j) = quad_handle.integrate(Functional::scalar_product(gradw_nTF, velocity_edge_basis->ancestor().function(j)), cell->edge(iTE));
                }
                for (size_t j = 0; j < n_local_velocity_cell_dofs; ++j)
                {
                    auto vT_on_F(Functional::trace(velocity_cell_basis->ancestor().function(j), edge_param));
                    gradw_nTF_vT(i, j) = quad_handle.integrate(Functional::scalar_product(gradw_nTF, vT_on_F), cell->edge(iTE));
                }
            }

            RT_RHS.topLeftCorner(n_local_velocity_highorder_dofs, n_local_velocity_cell_dofs) -= cell->edge_orientation(iTE) * velocity_highorder_basis->matrix() * gradw_nTF_vT * velocity_cell_basis->matrix().transpose();
            RT_RHS.block(0, local_offset, n_local_velocity_highorder_dofs, n_local_velocity_edge_dofs) = cell->edge_orientation(iTE) * velocity_highorder_basis->matrix() * gradw_nTF_vF * velocity_edge_basis->matrix().transpose();
        }

        Eigen::MatrixXd MFF = quad_handle.l2_product(*velocity_edge_basis, *velocity_edge_basis, cell->edge(iTE), sym);
        Eigen::MatrixXd MFHIGH = quad_handle.l2_product(*velocity_edge_basis, *velocity_highorder_basis, cell->edge(iTE), not_sym);

        for (size_t i = 0; i < n_local_velocity_edge_dofs; ++i)
        {
            for (size_t j = 0; j < n_local_velocity_edge_dofs; ++j)
            {
                triplets_M_BDRY_BDRY.emplace_back(m_stokes.local_offset_E(iT, iTE) + i, m_stokes.local_offset_E(iT, iTE) + j, MFF(i, j));
            }
            for (size_t l = 0; l < n_local_velocity_highorder_dofs; ++l)
            {
                triplets_M_BDRY_HIGH.emplace_back(m_stokes.local_offset_E(iT, iTE) + i, l, MFHIGH(i, l));
            }
        }
    }

    M_BDRY_BDRY.setFromTriplets(std::begin(triplets_M_BDRY_BDRY), std::end(triplets_M_BDRY_BDRY));
    M_BDRY_HIGH.setFromTriplets(std::begin(triplets_M_BDRY_HIGH), std::end(triplets_M_BDRY_HIGH));

    RT_loc = (ST + scalT * L_HIGH_HIGH).ldlt().solve(RT_RHS);
    Eigen::MatrixXd ATCONS = RT_loc.transpose() * ST * RT_loc;

    Eigen::MatrixXd DT = M_CELL_CELL.inverse() * M_CELL_HIGH * RT_loc;
    DT.topLeftCorner(n_local_velocity_cell_dofs, n_local_velocity_cell_dofs) -= Eigen::MatrixXd::Identity(n_local_velocity_cell_dofs, n_local_velocity_cell_dofs);

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_BDRY_BDRY);

    Eigen::MatrixXd M_BDRY_BDRY_INV_M_BDRY_HIGH = solver.solve(M_BDRY_HIGH);
    Eigen::MatrixXd DBDRY = M_BDRY_BDRY_INV_M_BDRY_HIGH * RT_loc;
    // Eigen::MatrixXd DBDRY = Eigen::MatrixXd::Zero(n_local_velocity_bdry_dofs, n_local_velocity_cell_dofs + n_local_velocity_bdry_dofs);
    // DBDRY.topLeftCorner(n_local_velocity_bdry_dofs, n_local_velocity_cell_dofs) = M_BDRY_BDRY_INV_M_BDRY_HIGH;
    DBDRY.bottomRightCorner(n_local_velocity_bdry_dofs, n_local_velocity_bdry_dofs) -= Eigen::MatrixXd::Identity(n_local_velocity_bdry_dofs, n_local_velocity_bdry_dofs);

    double hT = mesh_ptr->cell(iT)->diam();

    // Eigen::MatrixXd ATSTAB = pow(hT, -1) * DBDRY.transpose() * M_BDRY_BDRY * DBDRY;

    Eigen::MatrixXd ATSTAB = pow(hT, -2) * DT.transpose() * M_CELL_CELL * DT + pow(hT, -1) * DBDRY.transpose() * M_BDRY_BDRY * DBDRY;
    // Eigen::MatrixXd ATSTAB = DT.transpose() * quad_handle.h1_product(*velocity_cell_basis, *velocity_cell_basis, cell, sym) * DT + pow(hT, -1) * DBDRY.transpose() * M_BDRY_BDRY * DBDRY;

    AT_loc = Eigen::MatrixXd::Zero(n_local_velocity_bdry_dofs + n_local_velocity_cell_dofs + n_local_pressure_dofs, n_local_velocity_bdry_dofs + n_local_velocity_cell_dofs + n_local_pressure_dofs);

    AT_loc.topLeftCorner(n_local_velocity_bdry_dofs + n_local_velocity_cell_dofs, n_local_velocity_bdry_dofs + n_local_velocity_cell_dofs) = ATCONS + ATSTAB;
    AT_loc.bottomLeftCorner(n_local_pressure_dofs, n_local_velocity_bdry_dofs + n_local_velocity_cell_dofs) = -DT_RHS;
    AT_loc.topRightCorner(n_local_velocity_bdry_dofs + n_local_velocity_cell_dofs, n_local_pressure_dofs) = -DT_RHS.transpose();

    Eigen::MatrixXd pressure_mass_mat = quad_handle.l2_product(*pressure_basis, *pressure_basis, cell, sym);

    DivT[iT] = pressure_mass_mat.ldlt().solve(DT_RHS);
}

void Model::local_source_term(const size_t iT, Eigen::VectorXd &bT)
{
    bT = Eigen::VectorXd::Zero(m_stokes.local_cell_dofs(iT) + m_stokes.local_boundary_dofs(iT));

    auto quad_handle(m_stokes.get_quad_handle());
    bT.head(m_stokes.local_cell_dofs(iT)) = quad_handle.l2_product(m_src, *m_stokes.cell_basis(iT), mesh_ptr->cell(iT));
}

void Model::assemble(bool threading)
{
    // This method assembles the global matrices for the statically condensed system.
    // This static condensation procedure is developed in [Di-Pietro & Ern & Linke & Schieweck (2016)]

    std::cout << "\nAssembling the linear system\n";
    std::cout << "     (multi-threading = " << (threading ? "true" : "false") << ")\n";

    // Start the timer
    boost::timer::cpu_timer timer;
    timer.start();

    const size_t n_cells = mesh_ptr->n_cells();

    const size_t n_total_velocity_cell_dofs = m_stokes.total_cell_dofs();
    const size_t n_total_velocity_edge_dofs = m_stokes.total_edge_dofs();

    const size_t n_known_pressure_dofs = m_stokes.total_pressure_dofs() - n_cells; // only constant term on each cell unknown

    const size_t n_unknowns = n_total_velocity_edge_dofs + n_cells;

    // extra dof for lagrange multiplier to enforce zero integral condition

    inv_LTll_LTlg.resize(n_total_velocity_cell_dofs + n_known_pressure_dofs, n_unknowns + 1);
    LTgl.resize(n_unknowns + 1, n_total_velocity_cell_dofs + n_known_pressure_dofs);
    LTgg.resize(n_unknowns + 1, n_unknowns + 1);

    std::vector<Eigen::Triplet<double>> triplets_inv_LTll_LTlg;
    std::vector<Eigen::Triplet<double>> triplets_LTgl;
    std::vector<Eigen::Triplet<double>> triplets_LTgg;

    inv_LTll_rTl = Eigen::VectorXd::Zero(n_total_velocity_cell_dofs + n_known_pressure_dofs);
    std::vector<Eigen::MatrixXd> local_inv_LTll_LTlg(n_cells);
    std::vector<Eigen::MatrixXd> local_LTgl(n_cells);
    std::vector<Eigen::MatrixXd> local_LTgg(n_cells);

    rTg = Eigen::VectorXd::Zero(n_unknowns + 1);

    std::vector<Eigen::VectorXd> bT;

    RT.resize(n_cells);
    AT.resize(n_cells);
    bT.resize(n_cells);
    DivT.resize(n_cells);

    // Construct the local matrices using multithreading if threading is true
    std::function<void(size_t, size_t)> construct_all_local_contributions = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            // Eigen::VectorXd bT;
            local_stokes_operator(iT, AT[iT], RT[iT]);
            local_source_term(iT, bT[iT]);

            const size_t n_local_velocity_cell_dofs = m_stokes.local_cell_dofs(iT);
            const size_t n_local_velocity_bdry_dofs = m_stokes.local_boundary_dofs(iT);
            const size_t n_local_pressure_dofs = m_stokes.local_pressure_dofs(iT);

            const size_t n_local_known_dofs = n_local_velocity_cell_dofs + n_local_pressure_dofs - 1;
            const size_t n_local_velocity_dofs = n_local_velocity_cell_dofs + n_local_velocity_bdry_dofs;

            Eigen::MatrixXd local_LTll = Eigen::MatrixXd::Zero(n_local_known_dofs, n_local_known_dofs);
            local_LTll.topLeftCorner(n_local_velocity_cell_dofs, n_local_velocity_cell_dofs) = AT[iT].topLeftCorner(n_local_velocity_cell_dofs, n_local_velocity_cell_dofs);
            local_LTll.topRightCorner(n_local_velocity_cell_dofs, n_local_pressure_dofs - 1) = AT[iT].block(0, n_local_velocity_dofs + 1, n_local_velocity_cell_dofs, n_local_pressure_dofs - 1);
            local_LTll.bottomLeftCorner(n_local_pressure_dofs - 1, n_local_velocity_cell_dofs) = AT[iT].block(n_local_velocity_dofs + 1, 0, n_local_pressure_dofs - 1, n_local_velocity_cell_dofs);

            Eigen::MatrixXd local_LTlg = Eigen::MatrixXd::Zero(n_local_known_dofs, n_local_velocity_bdry_dofs + 1);
            local_LTlg.topLeftCorner(n_local_velocity_cell_dofs, n_local_velocity_bdry_dofs) = AT[iT].block(0, n_local_velocity_cell_dofs, n_local_velocity_cell_dofs, n_local_velocity_bdry_dofs);
            local_LTlg.bottomLeftCorner(n_local_pressure_dofs - 1, n_local_velocity_bdry_dofs) = AT[iT].block(n_local_velocity_dofs + 1, n_local_velocity_cell_dofs, n_local_pressure_dofs - 1, n_local_velocity_bdry_dofs);

            // local_LTgl[iT] = Eigen::MatrixXd::Zero(n_local_velocity_bdry_dofs + 1, n_local_known_dofs);
            // local_LTgl[iT].topLeftCorner(n_local_velocity_bdry_dofs, n_local_velocity_cell_dofs) = AT[iT].block(n_local_velocity_cell_dofs, 0, n_local_velocity_bdry_dofs, n_local_velocity_cell_dofs);
            // local_LTgl[iT].topRightCorner(n_local_velocity_bdry_dofs, n_local_pressure_dofs - 1) = AT[iT].block(n_local_velocity_cell_dofs, n_local_velocity_dofs + 1, n_local_velocity_bdry_dofs, n_local_pressure_dofs - 1);

            local_LTgl[iT] = local_LTlg.transpose(); // symmetry

            local_LTgg[iT] = AT[iT].block(n_local_velocity_cell_dofs, n_local_velocity_cell_dofs, n_local_velocity_bdry_dofs + 1, n_local_velocity_bdry_dofs + 1);

            Eigen::VectorXd local_RHS = Eigen::VectorXd::Zero(n_local_known_dofs);
            local_RHS.head(n_local_velocity_cell_dofs) = bT[iT].head(n_local_velocity_cell_dofs);

            Eigen::PartialPivLU<Eigen::MatrixXd> solver;
            solver.compute(local_LTll);

            Eigen::VectorXd local_inv_LTll_rTl = solver.solve(local_RHS);
            local_inv_LTll_LTlg[iT] = solver.solve(local_LTlg);

            inv_LTll_rTl.segment(m_stokes.global_offset_T(iT), n_local_velocity_cell_dofs) = local_inv_LTll_rTl.head(n_local_velocity_cell_dofs);
            inv_LTll_rTl.segment(n_total_velocity_cell_dofs + m_stokes.global_pressure_offset_T(iT) - iT, n_local_pressure_dofs - 1) = local_inv_LTll_rTl.tail(n_local_pressure_dofs - 1);
        }
    };

    // Running the local constructions in parallel
    parallel_for(n_cells, construct_all_local_contributions, threading);

    for (size_t iT = 0; iT < n_cells; ++iT)
    {
        Cell *cell = mesh_ptr->cell(iT);

        const size_t n_local_velocity_cell_dofs = m_stokes.local_cell_dofs(iT);
        const size_t n_local_velocity_bdry_dofs = m_stokes.local_boundary_dofs(iT);
        const size_t n_local_pressure_dofs = m_stokes.local_pressure_dofs(iT);

        for (size_t iTF = 0; iTF < cell->n_edges(); iTF++)
        {
            size_t iF = cell->edge(iTF)->global_index();
            rTg.segment(m_stokes.global_offset_E(iF), m_stokes.local_edge_dofs(iF)) += bT[iT].segment(n_local_velocity_cell_dofs + m_stokes.local_offset_E(iT, iTF), m_stokes.local_edge_dofs(iF));
            for (size_t i = 0; i < m_stokes.local_edge_dofs(iF); i++)
            {
                const size_t iGlobal = m_stokes.global_offset_E(iF) + i;
                const size_t iLocal = m_stokes.local_offset_E(iT, iTF) + i;
                for (size_t j = 0; j < n_local_velocity_cell_dofs; j++)
                {
                    const size_t jGlobal = m_stokes.global_offset_T(iT) + j;
                    triplets_inv_LTll_LTlg.emplace_back(jGlobal, iGlobal, local_inv_LTll_LTlg[iT](j, iLocal));
                    triplets_LTgl.emplace_back(iGlobal, jGlobal, local_LTgl[iT](iLocal, j));
                }
                for (size_t j = 0; j < n_local_pressure_dofs - 1; j++)
                {
                    size_t jGlobal = n_total_velocity_cell_dofs + m_stokes.global_pressure_offset_T(iT) - iT + j;
                    size_t jLocal = n_local_velocity_cell_dofs + j;
                    triplets_inv_LTll_LTlg.emplace_back(jGlobal, iGlobal, local_inv_LTll_LTlg[iT](jLocal, iLocal));
                    triplets_LTgl.emplace_back(iGlobal, jGlobal, local_LTgl[iT](iLocal, jLocal));
                }
                for (size_t jTF = 0; jTF < cell->n_edges(); jTF++)
                {
                    size_t jF = cell->edge(jTF)->global_index();
                    for (size_t j = 0; j < m_stokes.local_edge_dofs(jF); j++)
                    {
                        size_t jGlobal = m_stokes.global_offset_E(jF) + j;
                        const size_t jLocal = m_stokes.local_offset_E(iT, jTF) + j;
                        triplets_LTgg.emplace_back(iGlobal, jGlobal, local_LTgg[iT](iLocal, jLocal));
                    }
                }
                triplets_LTgg.emplace_back(iGlobal, n_unknowns - n_cells + iT, local_LTgg[iT](iLocal, n_local_velocity_bdry_dofs));
                triplets_LTgg.emplace_back(n_unknowns - n_cells + iT, iGlobal, local_LTgg[iT](n_local_velocity_bdry_dofs, iLocal));
            }
        }

        double pressure_integral = cell->measure() * m_stokes.pressure_basis(iT)->value(0, cell->center_mass());

        triplets_LTgg.emplace_back(n_total_velocity_edge_dofs + iT, n_unknowns, pressure_integral);
        triplets_LTgg.emplace_back(n_unknowns, n_total_velocity_edge_dofs + iT, pressure_integral);

        // triplets_LTgg.emplace_back(n_total_velocity_edge_dofs + iT, n_total_velocity_edge_dofs + iT, std::pow((cell->diam() / double(m_stokes.edge_degree() + 1)), m_stokes.edge_degree() + 2));
    }

    inv_LTll_LTlg.setFromTriplets(std::begin(triplets_inv_LTll_LTlg), std::end(triplets_inv_LTll_LTlg));
    LTgl.setFromTriplets(std::begin(triplets_LTgl), std::end(triplets_LTgl));
    LTgg.setFromTriplets(std::begin(triplets_LTgg), std::end(triplets_LTgg));

    // Get assembly time
    size_t assembly_time = timer.elapsed().wall;
    std::cout << "     Assembly time = " << double(assembly_time) * std::pow(10, -9) << "s\n";
}

Eigen::VectorXd Model::solve(const Eigen::VectorXd &UDir)
{
    std::cout << "\nSolving the linear system\n";

    // Start the timer
    boost::timer::cpu_timer timer;
    timer.start();

    const size_t n_cells = mesh_ptr->n_cells();

    const size_t n_total_velocity_cell_dofs = m_stokes.total_cell_dofs();
    const size_t n_total_velocity_edge_dofs = m_stokes.total_edge_dofs();
    const size_t n_total_pressure_dofs = m_stokes.total_pressure_dofs();

    size_t bdry_vector_dofs = UDir.size();

    Eigen::SparseMatrix<double> SysMat = (LTgg - LTgl * inv_LTll_LTlg);
    for (size_t ibF = 0; ibF < mesh_ptr->n_b_edges(); ++ibF)
    {
        const size_t iF = mesh_ptr->b_edge(ibF)->global_index();
        const size_t n_local_velocity_edge_dofs = m_stokes.local_edge_dofs(iF);
        const size_t offset = m_stokes.global_offset_E(iF);
        for (size_t i = 0; i < n_local_velocity_edge_dofs; ++i)
        {
            int index = int(offset + i);

            SysMat.row(index) *= 0.0;
        }
    }

    Eigen::VectorXd dirichlet_terms = SysMat.block(0, n_total_velocity_edge_dofs - bdry_vector_dofs, n_total_velocity_edge_dofs + n_cells + 1, bdry_vector_dofs) * UDir;

    for (size_t ibF = 0; ibF < mesh_ptr->n_b_edges(); ++ibF)
    {
        const size_t iF = mesh_ptr->b_edge(ibF)->global_index();
        const size_t n_local_velocity_edge_dofs = m_stokes.local_edge_dofs(iF);
        const size_t offset = m_stokes.global_offset_E(iF);
        for (size_t i = 0; i < n_local_velocity_edge_dofs; ++i)
        {
            int index = int(offset + i);

            SysMat.col(index) *= 0.0;
            SysMat.coeffRef(index, index) = 1.0;

            LTgl.row(index) *= 0.0;
        }
    }
    SysMat.prune(0.0);
    LTgl.prune(0.0);

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(SysMat);
    solver.factorize(SysMat);

    // Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;
    // solver.compute(SysMat);

    Eigen::VectorXd RHS = rTg - LTgl * inv_LTll_rTl - dirichlet_terms;

    Eigen::VectorXd condensed_terms = solver.solve(RHS);
    double residual_of_linear_system = (SysMat * condensed_terms - RHS).norm();
    std::cout << "     Residual of linear system = " << residual_of_linear_system << "\n";

    condensed_terms.segment(n_total_velocity_edge_dofs - bdry_vector_dofs, bdry_vector_dofs) = UDir;

    Eigen::VectorXd implicit_terms = inv_LTll_rTl - inv_LTll_LTlg * condensed_terms;

    Eigen::VectorXd U = Eigen::VectorXd::Zero(n_total_velocity_cell_dofs + n_total_velocity_edge_dofs + n_total_pressure_dofs);

    U.segment(0, n_total_velocity_cell_dofs) = implicit_terms.segment(0, n_total_velocity_cell_dofs);
    U.segment(n_total_velocity_cell_dofs, n_total_velocity_edge_dofs) = condensed_terms.segment(0, n_total_velocity_edge_dofs);

    for (size_t iT = 0; iT < n_cells; ++iT)
    {
        U(n_total_velocity_cell_dofs + n_total_velocity_edge_dofs + m_stokes.global_pressure_offset_T(iT)) = condensed_terms(n_total_velocity_edge_dofs + iT);

        U.segment(n_total_velocity_cell_dofs + n_total_velocity_edge_dofs + m_stokes.global_pressure_offset_T(iT) + 1, m_stokes.local_pressure_dofs(iT) - 1) = implicit_terms.segment(n_total_velocity_cell_dofs + m_stokes.global_pressure_offset_T(iT) - iT, m_stokes.local_pressure_dofs(iT) - 1);
    }

    // Get assembly time
    size_t solve_time = timer.elapsed().wall;
    std::cout << "     Solve time = " << double(solve_time) * std::pow(10, -9) << "s\n";

    return U;
}

void Model::plot(const Eigen::VectorXd &approx_Uvec)
{
    // std::ofstream vel_out("velocity_plot.csv");

    // // vel_out << "{";
    // for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++)
    // {
    //     Eigen::VectorXd rT_vec = RT[iT] * m_stokes.velocity_restr(approx_Uvec, iT);
    //     auto basis = m_stokes.highorder_basis(iT);

    //     std::function<Eigen::Vector2d(Eigen::Vector2d)> approx_sol = [&rT_vec, &basis](const Eigen::Vector2d &x) -> Eigen::Vector2d
    //     {
    //         Eigen::Vector2d value = Eigen::Vector2d::Zero();
    //         for (size_t i = 0; i < basis->dimension(); i++)
    //         {
    //             value += rT_vec(i) * basis->value(i, x);
    //         }
    //         return value;
    //     };

    //     QuadratureRule<Eigen::Vector2d> quadT = m_stokes.get_quad_handle().get_cell_quad(iT);

    //     for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
    //     {
    //         Eigen::Vector2d point(quadT[iqn].x);

    //         // if(iT !=  0)
    //         // {
    //         //     vel_out << ",";
    //         // }
    //         // Eigen::Vector2d point(mesh_ptr->cell(iT)->center_mass());
    //         // vel_out << "{{" << std::setprecision(12) << point(0) << "," << std::setprecision(12) << point(1)  << "},{" << approx_sol(point)(0) << "," << std::setprecision(12) << approx_sol(point)(1)  << "}}";

    //         vel_out << std::setprecision(12) << point.x() << "," << std::setprecision(12) << point.y() << "," << std::setprecision(12) << approx_sol(point).x() << "," << std::setprecision(12) << approx_sol(point).y() << "\n";

    //         // vel_out << std::setprecision(12) << std::setw(20) << std::left << point(0) << std::setprecision(12) << std::setw(20) << std::left << point(1) << std::setprecision(12) << std::setw(20) << std::left << approx_sol(point)(0) << std::setprecision(12) << std::setw(20) << std::left << approx_sol(point)(1) << std::endl;
    //     }
    // }
    // // vel_out << "}";
    // vel_out.close();

    std::vector<double> cell_enrich_values;
    for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++)
    {
        size_t L = m_stokes.cell_degree();
        int n_enrichments = m_stokes.local_cell_dofs(iT) - (L + 1) * (L + 2);
        cell_enrich_values.push_back(double(n_enrichments));
    }

    std::vector<Eigen::Vector2d> velocity_vertex_values;
    std::vector<double> pressure_vertex_values;

    for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++)
    {
        auto xV = mesh_ptr->vertex(iV)->coords();
        auto cList = mesh_ptr->vertex(iV)->get_cells();

        Eigen::Vector2d velocity_val = Eigen::Vector2d::Zero();
        double pressure_val = 0.0;
        for (size_t ilT = 0; ilT < cList.size(); ilT++)
        {
            size_t iT = cList[ilT]->global_index();

            Eigen::VectorXd local_potential = RT[iT] * m_stokes.velocity_restr(approx_Uvec, iT);

            for (size_t i = 0; i < m_stokes.local_highorder_dofs(iT); ++i)
            {
                velocity_val += local_potential(i) * m_stokes.highorder_basis(iT)->value(i, xV);
            }

            Eigen::VectorXd local_pressure = m_stokes.pressure_restr(approx_Uvec, iT);

            for (size_t i = 0; i < m_stokes.local_pressure_dofs(iT); ++i)
            {
                pressure_val += local_pressure(i) * m_stokes.pressure_basis(iT)->value(i, xV);
            }
        }
        velocity_val = velocity_val / (cList.size());
        pressure_val = pressure_val / (cList.size());

        velocity_vertex_values.push_back(velocity_val);
        pressure_vertex_values.push_back(pressure_val);
    }

    // vtk_writer plotdata(mesh_ptr, vertex_values);
    // plotdata.write_to_vtk("velocity.vtu");

    VtuWriter plotdata(mesh_ptr);
    plotdata.write_to_vtu("velocity.vtu", velocity_vertex_values);
    plotdata.write_to_vtu("pressure.vtu", pressure_vertex_values, false);
    plotdata.write_cells_to_vtu("cells2.vtu", cell_enrich_values);
}

std::vector<double> Model::compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::vector<Eigen::VectorXd> &elliptic_projectors, const std::vector<Eigen::VectorXd> &l2_velocity_projectors)
{
    std::cout << "\n[Scheme] Computing error\n";



    std::vector<double> pressure_l2_norm(mesh_ptr->n_cells());
    std::vector<double> velocity_h1_norm(mesh_ptr->n_cells());

    std::vector<double> pressure_diff_vec(mesh_ptr->n_cells());
    std::vector<double> pressure_exact_vec(mesh_ptr->n_cells());
    std::vector<double> energy_diff_vec(mesh_ptr->n_cells());
    std::vector<double> energy_exact_vec(mesh_ptr->n_cells());
    std::vector<double> L2_diff_vec(mesh_ptr->n_cells());
    std::vector<double> L2_exact_vec(mesh_ptr->n_cells());
    std::vector<double> H1_diff_vec(mesh_ptr->n_cells());
    std::vector<double> H1_exact_vec(mesh_ptr->n_cells());
    std::vector<double> divergence_vec(mesh_ptr->n_cells());
    // double discrete_pressure_avg = 0.0;
    // double interpolant_pressure_avg = 0.0;

    // for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
    // {
    //     Eigen::VectorXd local_interp_pressure = m_stokes.pressure_restr(interp_Uvec, iT);
    //     Eigen::VectorXd local_approx_pressure = m_stokes.pressure_restr(approx_Uvec, iT);

    //     auto quad_handle(m_stokes.get_quad_handle());
    //     QuadratureRule<Eigen::Vector2d> quadT(quad_handle.get_cell_quad(iT));

    //     for(int i = 0; i < local_interp_pressure.size(); ++i)
    //     {
    //         double basis_i_avg = 0.0;
    //         for(size_t iq = 0; iq < quadT.size(); ++iq)
    //         {
    //             basis_i_avg += quadT.weight(iq) * m_stokes.pressure_basis(iT)->value(i, quadT.point(iq));
    //         }

    //         discrete_pressure_avg += local_approx_pressure[i] * basis_i_avg;
    //         interpolant_pressure_avg += local_interp_pressure[i] * basis_i_avg;
    //     }
    // }
    // std::cout << "     Discrete pressure avg = " << discrete_pressure_avg << "\n";
    // std::cout << "     Interpolant pressure avg = " << interpolant_pressure_avg << "\n\n";

    std::function<void(size_t, size_t)> compute_all_errors = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; ++iT)
        {
            Eigen::VectorXd local_interp_velocity = m_stokes.velocity_restr(interp_Uvec, iT);
            Eigen::VectorXd local_interp_pressure = m_stokes.pressure_restr(interp_Uvec, iT);

            Eigen::VectorXd local_approx_velocity = m_stokes.velocity_restr(approx_Uvec, iT);
            Eigen::VectorXd local_approx_pressure = m_stokes.pressure_restr(approx_Uvec, iT);

            Eigen::VectorXd local_difference_velocity = local_interp_velocity - local_approx_velocity;
            Eigen::VectorXd local_difference_pressure = local_interp_pressure - local_approx_pressure;

            const size_t local_velocity_dofs = m_stokes.local_cell_dofs(iT) + m_stokes.local_boundary_dofs(iT);

            energy_diff_vec[iT] = local_difference_velocity.transpose() * AT[iT].topLeftCorner(local_velocity_dofs, local_velocity_dofs) * local_difference_velocity;
            energy_exact_vec[iT] = local_interp_velocity.transpose() * AT[iT].topLeftCorner(local_velocity_dofs, local_velocity_dofs) * local_interp_velocity;

            auto quad_handle(m_stokes.get_quad_handle());

            Eigen::MatrixXd ST = quad_handle.h1_product(*m_stokes.highorder_basis(iT), *m_stokes.highorder_basis(iT), mesh_ptr->cell(iT), true);
            Eigen::MatrixXd MT = quad_handle.l2_product(*m_stokes.highorder_basis(iT), *m_stokes.highorder_basis(iT), mesh_ptr->cell(iT), true);

            Eigen::MatrixXd MPressure = quad_handle.l2_product(*m_stokes.pressure_basis(iT), *m_stokes.pressure_basis(iT), mesh_ptr->cell(iT), true);

            Eigen::VectorXd rT_minus_sigma_T = RT[iT] * local_approx_velocity - elliptic_projectors[iT];
            Eigen::VectorXd sigma_T = elliptic_projectors[iT];

            Eigen::VectorXd rT_minus_pi_T = RT[iT] * local_approx_velocity - l2_velocity_projectors[iT];
            Eigen::VectorXd pi_T = l2_velocity_projectors[iT];

            H1_diff_vec[iT] = rT_minus_sigma_T.transpose() * ST * rT_minus_sigma_T;
            H1_exact_vec[iT] = sigma_T.transpose() * ST * sigma_T;

            velocity_h1_norm[iT] = (RT[iT] * local_approx_velocity).transpose() * ST * RT[iT] * local_approx_velocity;
            pressure_l2_norm[iT] = local_approx_pressure.transpose() * MPressure * local_approx_pressure;

            L2_diff_vec[iT] = rT_minus_pi_T.transpose() * MT * rT_minus_pi_T;
            L2_exact_vec[iT] = pi_T.transpose() * MT * pi_T;

            pressure_diff_vec[iT] = local_difference_pressure.transpose() * MPressure * local_difference_pressure;
            pressure_exact_vec[iT] = local_interp_pressure.transpose() * MPressure * local_interp_pressure;

            divergence_vec[iT] = (DivT[iT] * local_approx_velocity).transpose() * MPressure * (DivT[iT] * local_approx_velocity);
        }
    };
    parallel_for(mesh_ptr->n_cells(), compute_all_errors, true);

    std::cout << std::setprecision(16) << "     Velocity H1 norm = " << std::sqrt(std::accumulate(velocity_h1_norm.begin(), velocity_h1_norm.end(), 0.0)) << "\n";
    std::cout << std::setprecision(16) << "     Pressure L2 norm = " << std::sqrt(std::accumulate(pressure_l2_norm.begin(), pressure_l2_norm.end(), 0.0)) << "\n\n";

    // double pressure_diff = std::accumulate(pressure_diff_vec.begin(), pressure_diff_vec.end(), 0.0);
    // double pressure_exact = std::accumulate(pressure_exact_vec.begin(), pressure_exact_vec.end(), 0.0);

    double pressure_diff = std::abs(std::sqrt(std::accumulate(pressure_l2_norm.begin(), pressure_l2_norm.end(), 0.0)) - 15.61002388089276);

    double energy_diff = std::accumulate(energy_diff_vec.begin(), energy_diff_vec.end(), 0.0);
    double energy_exact = std::accumulate(energy_exact_vec.begin(), energy_exact_vec.end(), 0.0);

    double L2_diff = std::accumulate(L2_diff_vec.begin(), L2_diff_vec.end(), 0.0);
    double L2_exact = std::accumulate(L2_exact_vec.begin(), L2_exact_vec.end(), 0.0);

    // double H1_diff = std::accumulate(H1_diff_vec.begin(), H1_diff_vec.end(), 0.0);
    // double H1_exact = std::accumulate(H1_exact_vec.begin(), H1_exact_vec.end(), 0.0);

    double H1_diff = std::abs(std::sqrt(std::accumulate(velocity_h1_norm.begin(), velocity_h1_norm.end(), 0.0)) - 5.777266571568622);

    double divergence = std::accumulate(divergence_vec.begin(), divergence_vec.end(), 0.0);

    // return {std::sqrt(L2_diff / L2_exact), std::sqrt(H1_diff / H1_exact), std::sqrt(energy_diff / energy_exact), std::sqrt(pressure_diff / pressure_exact), std::sqrt(divergence)};
    // return {std::sqrt(L2_diff / L2_exact), std::sqrt(H1_diff / H1_exact), std::sqrt(energy_diff / energy_exact), std::sqrt(pressure_diff)};
    return {std::sqrt(L2_diff), std::sqrt(H1_diff), std::sqrt(energy_diff), std::sqrt(pressure_diff), std::sqrt(divergence)};
}