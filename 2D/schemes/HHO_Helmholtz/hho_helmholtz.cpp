#include "hho_helmholtz.hpp"

// internal libraries
#include "parallel_for.hpp"
#include "vtu_writer.hpp"
#include "QuadratureRule.hpp"
#include "basis.hpp"
#include "function.hpp"

// mesh intersection
#include "../IntersectMesh/IntersectMesh.hpp"
#include "TestCase.hpp"

// boost libraries
#include <boost/program_options.hpp> // boost::program_options
#include <boost/timer/timer.hpp>     // boost::cpu_timer

// std libraries
#include <iostream> // std::cout

using namespace PolyMesh2D::HHOHELMHOLTZ;

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

    TestCase TC(params.test_case, params.bdry_case);

    PolyMesh2D::MeshCutter mesh_cutter(straight_mesh.get(), TC.get_level_set(), TC.get_boundary_param());
    std::unique_ptr<Mesh> curved_mesh = mesh_cutter.cut_mesh();

    straight_mesh.reset(); // delete straight mesh

    // Reordering
    reorder_edges(curved_mesh.get());

    std::ofstream mesh_out("mesh_plot.dat");
    curved_mesh->plot_mesh(&mesh_out, 20);
    mesh_out.close();
    // exit(0);

    assert(curved_mesh->test());

    // Print the data
    std::cout << "[Scheme] Data:\n";
    std::cout << "     No. cells = " << curved_mesh->n_cells() << ", No. edges = " << curved_mesh->n_edges() << ", No. vertices = " << curved_mesh->n_vertices() << "\n";
    std::cout << "     TestCase = " << params.test_case << "\n";
    std::cout << "     Mesh = " << params.mesh_name << "\n";
    std::cout << "     Degrees: edge = " << params.edge_degree << "; cell = " << params.cell_degree << "\n";
    std::cout << "     Using threads = " << (params.use_threads ? "true" : "false") << std::endl;

    HybridCore hho(curved_mesh.get(), params.cell_degree, params.edge_degree, params.use_threads, params.orthonormalise);

    ScalarFunction2D source(TC.src());
    ScalarFunction2D exact_sol(TC.sol());

    Model model(hho, source);

    model.assemble(params.use_threads);
    Eigen::VectorXd approx_Uvec(model.solve());
    Eigen::VectorXd interp_Uvec(hho.interpolate(exact_sol.get_value()));

    auto errors = model.compute_errors(approx_Uvec, interp_Uvec, exact_sol);

    double l2_error = errors[0];
    double h1_error = errors[1];
    double energy_error = errors[2];

    std::cout << "     L2 Error = " << l2_error << "\n";
    std::cout << "     H1 Error = " << h1_error << "\n";
    std::cout << "     Energy Error = " << energy_error << "\n";

    model.plot(approx_Uvec, interp_Uvec, params.plot_file, exact_sol);

    std::cout << "\n[Scheme] Writing solution to file\n";

    std::ofstream out("results.txt");
    out << "EdgeDegree: " << params.edge_degree << "\n";
    out << "CellDegree: " << params.cell_degree << "\n";
    out << "L2Error: " << l2_error << "\n";
    out << "H1Error: " << h1_error << "\n";
    out << "EnergyError: " << energy_error << "\n";
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
    desc.add_options()("help,h", "Produce help message")("mesh,m", po::value<std::string>(), "Set the mesh")("testcase,t", po::value<unsigned>(), "Set the exact solution")("celldegree,l", po::value<unsigned>(), "Set the degree of the cell polynomials")("edgedegree,k", po::value<unsigned>(), "Set the degree of the edge polynomials")("plot,p", po::value<std::string>(), "Plot to file")("use_threads,u", po::value<bool>(), "Using multithreading")("orthonormalise,o", po::value<bool>(), "Orthonormalise the basis functions")("boundary,b", po::value<char>(), "Set the boundary ('E' = ellipse, 'C' = circle)");

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
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "../../../typ2_meshes/mesh2_2_transformed.typ2");

    // Get test case
    test_case = (vm.count("testcase") ? vm["testcase"].as<unsigned>() : 1);

    if (test_case != 1 && test_case != 2 && test_case != 3)
    {
        std::cout << "Invalid choice of test case\n";
        exit(1);
    }

    // Get test case
    bdry_case = (vm.count("boundary") ? vm["boundary"].as<char>() : 'C');

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

Model::Model(const HybridCore &hho, const ScalarFunction2D &src) : m_hho(hho), m_src(src), mesh_ptr(m_hho.get_mesh()) {}

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

            local_helmholtz_operator(iT, AT[iT], PT[iT]);
            local_source_term(iT, bT[iT]);

            // Get each contribution of AT
            Eigen::MatrixXd ATT = AT[iT].topLeftCorner(n_local_cell_dofs, n_local_cell_dofs);
            Eigen::MatrixXd ATF = AT[iT].topRightCorner(n_local_cell_dofs, n_local_bdry_dofs);
            // Eigen::MatrixXd AFT = AT[iT].bottomLeftCorner(n_bdry_dofs, n_local_cell_dofs);
            Eigen::MatrixXd AFT = ATF.transpose(); // symmetry
            Eigen::MatrixXd AFF = AT[iT].bottomRightCorner(n_local_bdry_dofs, n_local_bdry_dofs);

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

Eigen::VectorXd Model::solve()
{
    std::cout << "\n[Scheme] Solving the linear system\n";

    // Start the timer
    boost::timer::cpu_timer assembly_timer;
    assembly_timer.start();

    size_t n_fixed_dofs = 0;

    for (size_t iF = 0; iF < mesh_ptr->n_edges(); ++iF)
    {
        if (!(mesh_ptr->edge(iF)->is_boundary()))
        {
            continue;
        }
        n_fixed_dofs += m_hho.local_edge_dofs(iF);
    }
    const size_t n_unknowns = m_hho.total_edge_dofs() - n_fixed_dofs;

    Eigen::SparseMatrix<double> SysMat = GlobMat.topLeftCorner(n_unknowns, n_unknowns);
    Eigen::VectorXd UDir = Eigen::VectorXd::Zero(n_fixed_dofs);

    Eigen::VectorXd B = GlobRHS.segment(0, n_unknowns) - GlobMat.topRightCorner(n_unknowns, n_fixed_dofs) * UDir;

    // Solve the statically condesed system using BiCGSTAB
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(SysMat);

    // xF = INV( AII - AIT * INV_ATT * ATI ) * ( BI - AIT * INV_ATT * BT - (AID - AIT * INV_ATT * ATD * UD) )
    Eigen::VectorXd xF = solver.solve(B);

    // Print solver iterations and estimated error
    std::cout << "     [solver] #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;
    // _iterations = solver.iterations();

    // Find the residual of the system
    // solving_error = (SysMat * xF - B).norm();

    // Set each part of the solution
    Eigen::VectorXd Xh = Eigen::VectorXd::Zero(m_hho.total_cell_dofs() + m_hho.total_edge_dofs());
    Xh.tail(n_fixed_dofs) = UDir;
    Xh.segment(m_hho.total_cell_dofs(), n_unknowns) = xF;

    //                         = INV_ATT_BT - INV_ATT * ATF * UF
    Xh.head(m_hho.total_cell_dofs()) = ScRHS - ScBeMat * Xh.tail(m_hho.total_edge_dofs());

    std::cout << "     Solving time = " << assembly_timer.elapsed().wall * pow(10, -9) << "s\n";

    return Xh;
}

std::vector<double> Model::compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const ScalarFunction2D &sol)
{
    std::cout << "\n[Scheme] Computing error\n";

    double energy_diff = 0.0;
    double energy_exact = 0.0;

    double L2_diff = 0.0;
    double L2_exact = 0.0;

    double H1_diff = 0.0;
    double H1_exact = 0.0;

    for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
    {
        Eigen::VectorXd ITk = m_hho.restr(interp_Uvec, iT);
        Eigen::VectorXd diff_vec = m_hho.restr(approx_Uvec, iT) - ITk;

        energy_diff += diff_vec.transpose() * AT[iT] * diff_vec;
        energy_exact += ITk.transpose() * AT[iT] * ITk;

        Eigen::VectorXd pT_vec = PT[iT] * m_hho.restr(approx_Uvec, iT);
        auto basis = m_hho.highorder_basis(iT);
        QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

        std::function<PolyMesh2D::Functional::RowVector(PolyMesh2D::Functional::ColVector)> approx_grad = [&pT_vec, &basis](const PolyMesh2D::Functional::ColVector &x) -> PolyMesh2D::Functional::RowVector
        {
            PolyMesh2D::Functional::RowVector value = PolyMesh2D::Functional::RowVector::Zero();
            for (size_t i = 0; i < basis->dimension(); i++)
            {
                value += pT_vec(i) * basis->derivative(i, x);
            }
            return value;
        };

        for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
        {
            H1_diff += quadT[iqn].w * (approx_grad(quadT[iqn].x) - sol.derivative(quadT[iqn].x)).squaredNorm();
            H1_exact += quadT[iqn].w * sol.derivative(quadT[iqn].x).squaredNorm();
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

        for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
        {
            L2_diff += quadT[iqn].w * std::pow(approx_sol(quadT[iqn].x) - sol.value(quadT[iqn].x), 2);
            L2_exact += quadT[iqn].w * std::pow(sol.value(quadT[iqn].x), 2);
        }
    }

    return {std::sqrt(L2_diff / L2_exact), std::sqrt(H1_diff / H1_exact), std::sqrt(energy_diff / energy_exact)};
}

void Model::plot(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::string &plot_file, const ScalarFunction2D &sol)
{
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
    }
}

void Model::local_helmholtz_operator(const size_t iT, Eigen::MatrixXd &AT, Eigen::MatrixXd &PT)
{
    const Cell *cell = mesh_ptr->cell(iT);
    const size_t n_local_edges = cell->n_edges();
    const size_t n_local_cell_dofs = m_hho.local_cell_dofs(iT);
    const size_t n_local_highorder_dofs = m_hho.local_highorder_dofs(iT);
    const size_t n_local_bdry_dofs = m_hho.local_boundary_dofs(iT);
    const size_t n_local_dofs = n_local_cell_dofs + n_local_bdry_dofs;

    Eigen::MatrixXd M_CELL_CELL = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_cell_dofs);
    Eigen::MatrixXd M_CELL_HIGH = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_highorder_dofs);

    QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);

    for (size_t i = 0; i < n_local_cell_dofs; ++i)
    {
        for (size_t j = 0; j < n_local_highorder_dofs; ++j)
        {
            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                M_CELL_HIGH(i, j) += quadT[iqn].w * m_hho.cell_basis(iT)->value(i, quadT[iqn].x) * m_hho.highorder_basis(iT)->value(j, quadT[iqn].x);
                if(j < n_local_cell_dofs)
                {
                    M_CELL_CELL(i, j) += quadT[iqn].w * m_hho.cell_basis(iT)->value(i, quadT[iqn].x) * m_hho.cell_basis(iT)->value(j, quadT[iqn].x);
                }
            }
            // if ((i < j) && (j < n_local_cell_dofs))
            // {
            //     MTT(j, i) = MTT(i, j);
            // }
        }
    }

    // Calculate stiffness matrix
    Eigen::MatrixXd ST_high_high = Eigen::MatrixXd::Zero(n_local_highorder_dofs, n_local_highorder_dofs);
    Eigen::MatrixXd ST_high_cell = Eigen::MatrixXd::Zero(n_local_highorder_dofs, n_local_cell_dofs);
    Eigen::MatrixXd ST_cell_cell = Eigen::MatrixXd::Zero(n_local_cell_dofs, n_local_cell_dofs);

    for (size_t i = 0; i < n_local_highorder_dofs; ++i)
    {
        for (size_t j = 0; j < n_local_highorder_dofs; ++j)
        {
            for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
            {
                ST_high_high(i, j) += quadT[iqn].w * m_hho.highorder_basis(iT)->derivative(i, quadT[iqn].x) * (m_hho.highorder_basis(iT)->derivative(j, quadT[iqn].x)).transpose();
                if(j < n_local_cell_dofs)
                {
                    ST_high_cell(i, j) += quadT[iqn].w * m_hho.highorder_basis(iT)->derivative(i, quadT[iqn].x) * (m_hho.cell_basis(iT)->derivative(j, quadT[iqn].x)).transpose();
                    if(i < n_local_cell_dofs)
                    {
                        ST_cell_cell(i, j) += quadT[iqn].w * m_hho.cell_basis(iT)->derivative(i, quadT[iqn].x) * (m_hho.cell_basis(iT)->derivative(j, quadT[iqn].x)).transpose();
                    }
                }
            }
            // ST(j, i) = ST(i, j);
        }
    }

    // Get L_T matrix to impose average condition on P_T later
    Eigen::MatrixXd LTtLT = M_CELL_HIGH.col(0) * M_CELL_HIGH.row(0);
    double scalT = ST.norm() / LTtLT.norm();

    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_BDRY;
    std::vector<Eigen::Triplet<double>> triplets_M_BDRY_HIGH;

    Eigen::SparseMatrix<double> M_BDRY_BDRY(n_local_bdry_dofs, n_local_bdry_dofs);
    Eigen::SparseMatrix<double> M_BDRY_HIGH(n_local_bdry_dofs, n_local_highorder_dofs);

    // Set cell - cell term of RHS of P_T
    Eigen::MatrixXd BP = Eigen::MatrixXd::Zero(n_local_highorder_dofs, n_local_dofs);
    BP.topLeftCorner(n_local_highorder_dofs, n_local_cell_dofs) = ST_high_cell + scalT * LTtLT;

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
                for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
                {
                    MFF_i_j += quadF[iqn].w * m_hho.edge_basis(iF)->value(i, quadF[iqn].x) * m_hho.edge_basis(iF)->value(j, quadF[iqn].x);
                }
                triplets_M_BDRY_BDRY.emplace_back(offset_F + i, offset_F + j, MFF_i_j);
                // if (i < j)
                // {
                //     triplets_M_BDRY_BDRY.emplace_back(offset_F + j, offset_F + i, MFF_i_j);
                // }
            }
            for (size_t l = 0; l < n_local_highorder_dofs; ++l)
            {
                double MFHIGH_i_l = 0.0;
                for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
                {
                    MFHIGH_i_l += quadF[iqn].w * m_hho.edge_basis(iF)->value(i, quadF[iqn].x) * m_hho.highorder_basis(iT)->value(l, edge_param.value(quadF[iqn].x));
                }
                triplets_M_BDRY_HIGH.emplace_back(offset_F + i, l, MFHIGH_i_l);
            }
        }

        for (size_t iqn = 0; iqn < quadF.size(); ++iqn)
        {
            for (size_t i = 0; i < n_local_highorder_dofs; ++i)
            {
                double grad_w_dot_nTF = m_hho.highorder_basis(iT)->derivative(i, edge_param.value(quadF[iqn].x)) * (cell->edge_normal(iTF, quadF[iqn].x));
                for (size_t j = 0; j < n_local_cell_dofs; ++j)
                {
                    BP(i, j) -= quadF[iqn].w * grad_w_dot_nTF * m_hho.cell_basis(iT)->value(j, edge_param.value(quadF[iqn].x));
                }
                for (size_t l = 0; l < n_local_edge_dofs; ++l)
                {
                    BP(i, offset_TF + l) += quadF[iqn].w * grad_w_dot_nTF * m_hho.edge_basis(iF)->value(l, quadF[iqn].x);
                }
            }
        }
    }

    M_BDRY_BDRY.setFromTriplets(std::begin(triplets_M_BDRY_BDRY), std::end(triplets_M_BDRY_BDRY));
    M_BDRY_HIGH.setFromTriplets(std::begin(triplets_M_BDRY_HIGH), std::end(triplets_M_BDRY_HIGH));

    PT = (ST_high_high + scalT * LTtLT).ldlt().solve(BP);
    Eigen::MatrixXd ATCONS = PT.transpose() * ST_high_high * PT;

    Eigen::MatrixXd DT = M_CELL_CELL.inverse() * M_CELL_HIGH * PT;
    DT.topLeftCorner(n_local_cell_dofs, n_local_cell_dofs) -= Eigen::MatrixXd::Identity(n_local_cell_dofs, n_local_cell_dofs);

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(M_BDRY_BDRY);

    Eigen::MatrixXd M_BDRY_BDRY_INV_M_BDRY_HIGH = solver.solve(M_BDRY_HIGH);
    Eigen::MatrixXd STAB_OPERATOR = M_BDRY_BDRY_INV_M_BDRY_HIGH * PT;
    STAB_OPERATOR.bottomRightCorner(n_local_bdry_dofs, n_local_bdry_dofs) -= Eigen::MatrixXd::Identity(n_local_bdry_dofs, n_local_bdry_dofs);

    Eigen::MatrixXd ATSTAB = DT.transpose() * ST_cell_cell * DT + (1.0 / cell->diam()) * STAB_OPERATOR.transpose() * M_BDRY_BDRY * STAB_OPERATOR;

    AT = ATCONS + ATSTAB;
}

void Model::local_source_term(const size_t iT, Eigen::VectorXd &bT)
{
    bT = Eigen::VectorXd::Zero(m_hho.local_cell_dofs(iT) + m_hho.local_boundary_dofs(iT));

    QuadratureRule<Eigen::Vector2d> quadT = m_hho.quadT(iT);
    for (size_t iqn = 0; iqn < quadT.size(); ++iqn)
    {
        double source_weight = quadT[iqn].w * m_src.value(quadT[iqn].x);
        for (size_t i = 0; i < m_hho.local_cell_dofs(iT); ++i)
        {
            bT(i) += source_weight * m_hho.cell_basis(iT)->value(i, quadT[iqn].x);
        }
    }
}
