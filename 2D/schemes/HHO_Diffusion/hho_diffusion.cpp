#include "hho_diffusion.hpp"

// internal libraries
#include "parallel_for.hpp"
#include "vtu_writer.hpp"
#include "QuadratureRule.hpp"
#include "basis.hpp"
#include "function.hpp"

#include "../HHO_Poisson/TestCase.hpp"

// mesh intersection
#include "../IntersectCurvedMesh/IntersectCurvedMesh.hpp"
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
    std::unique_ptr<Mesh> curved_mesh = mesh_cutter_outer.cut_mesh();

    straight_mesh.reset(); // delete straight mesh

    double R = 0.6;

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
    PolyMesh2D::MeshIntersect mesh_cutter(curved_mesh.get(), level_set, bdry_param);
    mesh_cutter.cut_mesh();

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

    double beta1 = 1E-6;
    double beta2 = 1.0;

    double a = 2.5;

    std::function<Eigen::Matrix2d(PolyMesh2D::CurvedMesh::Cell *)> diffusion = [&beta1, &beta2, &R](PolyMesh2D::CurvedMesh::Cell *cell) -> Eigen::Matrix2d
    {
        Eigen::Matrix2d K = Eigen::Matrix2d::Identity(2, 2);

        if (cell->center_mass().norm() < R)
        {
            // r < R
            K(0, 1) = 1.0 - beta1;
            K(1, 0) = 1.0 - beta1;
        }
        else
        {
            // r > R
            K(0, 1) = 1.0 - beta2;
            K(1, 0) = 1.0 - beta2;
        }

        return K;
    };

    std::function<double(Eigen::Vector2d, PolyMesh2D::CurvedMesh::Cell *)> src = [&a, &R](const Eigen::Vector2d &x, PolyMesh2D::CurvedMesh::Cell *cell) -> double
    {
        return 1.0;
    };

    Model model(hho, src, diffusion);

    size_t n_fixed_dofs = 0;
    for (auto &b_edge : curved_mesh->get_b_edges())
    {
        n_fixed_dofs += hho.local_edge_dofs(b_edge->global_index());
    }

    model.assemble(params.use_threads);

    Eigen::VectorXd UDir(Eigen::VectorXd::Zero(n_fixed_dofs));
    Eigen::VectorXd approx_Uvec(model.solve(UDir));

    auto errors = model.compute_errors(approx_Uvec);

    double integral = errors[0];
    double integral_error = errors[1];
    double h1_norm = errors[2];
    double h1_norm_error = errors[3];

    model.plot(approx_Uvec, params.plot_file);

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
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "/home/liam/github/codes/liamyemm/PolyMesh/2D/typ2_meshes/mesh2_2_transformed.typ2");

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

            local_diffusion_operator(iT);
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

std::vector<double> Model::compute_errors(const Eigen::VectorXd &approx_Uvec)
{
    std::vector<double> pT_average_vec(mesh_ptr->n_cells());
    std::vector<double> pT_H1_norm_vec(mesh_ptr->n_cells());

    std::function<void(size_t, size_t)> compute_average = [&](size_t start, size_t end) -> void
    {
    for (size_t iT = start; iT < end; iT++)
    {
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
    }
    };
    parallel_for(mesh_ptr->n_cells(), compute_average, true);
    double pT_average = 0.0;
    double pT_H1_norm = 0.0;

    for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
    {
        pT_average += pT_average_vec[iT];
        pT_H1_norm += pT_H1_norm_vec[iT];
    }
    // return {std::abs(pT_average - 0.47342316790031907514), std::abs(flux + Math::PI)};
    // return {pT_average, std::abs(pT_average - 0.46006947132345432649), std::sqrt(pT_H1_norm), std::abs(std::sqrt(pT_H1_norm) - 0.80699765696519920599)};
    return {pT_average, std::abs(pT_average - 0.46006947132349923502), std::sqrt(pT_H1_norm), std::abs(std::sqrt(pT_H1_norm) - 0.8069976569614020212)};
}

void Model::plot(const Eigen::VectorXd &approx_Uvec, const std::string &plot_file)
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
        Eigen::VectorXd potential = Eigen::VectorXd::Zero(mesh_ptr->n_vertices());

        for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++)
        {
            auto xV = mesh_ptr->vertex(iV)->coords();
            auto cList = mesh_ptr->vertex(iV)->get_cells();
            for (size_t ilT = 0; ilT < cList.size(); ilT++)
            {
                size_t iT = cList[ilT]->global_index();

                Eigen::VectorXd local_potential = PT[iT] * m_hho.restr(approx_Uvec, iT);

                for (size_t i = 0; i < m_hho.local_highorder_dofs(iT); i++)
                {
                    potential(iV) += local_potential(i) * m_hho.highorder_basis(iT)->value(i, xV);
                }
            }
            potential(iV) = potential(iV) / (cList.size());
        }

        VtuWriter plotdata(mesh_ptr);
        plotdata.write_to_vtu("potential-" + plot_file + ".vtu", potential);
    }
}

void Model::local_diffusion_operator(const size_t iT)
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