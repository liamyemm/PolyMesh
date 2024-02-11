#include "xvem_poisson.h"

// internal libraries
#include "parallel_for.hpp"
#include "vtu_writer.hpp"
#include "QuadratureRule.hpp"
#include "basis.hpp"
#include "function.hpp"

#include "MeshBuilder2D.hpp"
#include "MeshReaderTyp2.hpp"

#include "mesh_convert.hpp"
#include "PoissonSingularity.h"

// boost libraries
#include <boost/program_options.hpp> // boost::program_options
#include <boost/timer/timer.hpp>     // boost::cpu_timer

// #include "affine_shift.hpp"

#include <numeric> // std::accumulate

// std libraries
#include <iostream> // std::cout

using namespace PolyMesh2D::XVEM;

bool locally_enrich(PolyMesh2D::CurvedMesh::Cell * cell, double tol)
{
    bool enrich = false;
    for(auto & v : cell->get_vertices())
    {
        enrich = enrich || v->coords().norm() < tol;
    }
    return enrich;
}

bool locally_enrich(PolyMesh2D::CurvedMesh::Edge * edge, double tol)
{
    bool enrich = false;
    for(auto & c : edge->get_cells())
    {
        enrich = enrich || locally_enrich(c, tol);
    }
    return enrich;
}

int main(const int argc, const char **argv)
{
    ModelParameters params(argc, argv);

    // PolyMesh2D::affine_shift("../../../typ2_meshes/mesh2_0.typ2", "../../../typ2_meshes/mesh2_0_transformed.typ2", 2.0 * Eigen::Matrix2d::Identity(), Eigen::Vector2d(-1.0, -1.0));

    // exit(1);

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

    std::unique_ptr<Mesh> curved_mesh = PolyMesh2D::MeshTransform::mesh_convert(straight_mesh.get());

    straight_mesh.reset(); // delete straight mesh

    // Reordering
    // reorder_edges(curved_mesh.get());
    // exit(0);

    assert(curved_mesh->test());

    // for (auto &e : curved_mesh->get_edges())
    // {
    //     double x = e->center_mass()(0);
    //     double y = e->center_mass()(1);
    //     if (std::abs(y) < 1E-14 && x > 0)
    //     {
    //         e->set_boundary(true);
    //         for (auto &v : e->get_vertices())
    //         {
    //             v->set_boundary(true);
    //         }
    //         for (auto &c : e->get_cells())
    //         {
    //             c->set_boundary(true);
    //         }
    //     }
    // }

    // curved_mesh->reset_boundary(); // empty all boundary and internal vectors, so we can reset them without duplicating.

    // for (auto &cell : curved_mesh->get_cells())
    // {
    //     if (cell->is_boundary())
    //     {
    //         curved_mesh->add_b_cell(cell);
    //     }
    //     else
    //     {
    //         curved_mesh->add_i_cell(cell);
    //     }
    // }

    // for (auto &edge : curved_mesh->get_edges())
    // {
    //     if (edge->is_boundary())
    //     {
    //         curved_mesh->add_b_edge(edge);
    //     }
    //     else
    //     {
    //         curved_mesh->add_i_edge(edge);
    //     }
    // }

    // for (auto &vert : curved_mesh->get_vertices())
    // {
    //     if (vert->is_boundary())
    //     {
    //         curved_mesh->add_b_vertex(vert);
    //     }
    //     else
    //     {
    //         curved_mesh->add_i_vertex(vert);
    //     }
    // }

    std::ofstream mesh_out("mesh_plot.dat");
    curved_mesh->plot_mesh(&mesh_out, 20);
    mesh_out.close();

    // Print the data
    std::cout << "[Scheme] Data:\n";
    std::cout << "     No. cells = " << curved_mesh->n_cells() << ", No. edges = " << curved_mesh->n_edges() << ", No. vertices = " << curved_mesh->n_vertices() << "\n";
    // std::cout << "     TestCase = " << params.test_case << "\n";
    // std::cout << "     Mesh = " << params.mesh_name << "\n";
    std::cout << "     Degrees: k = " << params.edge_degree << "; l = " << params.cell_degree << "\n";
    std::cout << "     Using threads = " << (params.use_threads ? "true" : "false") << std::endl;

    PolyMesh2D::XVEMCore vem(curved_mesh.get(), params.cell_degree, params.edge_degree, params.use_threads, params.orthonormalise);

    PolyMesh2D::Functional::Function<2, 1> singularity = PolyMesh2D::PoissonSingularity(2.0 / 3.0, -Math::PI / 2.0, Eigen::Vector2d(0.0, 0.0));

    // PolyMesh2D::Functional::Function<2, 1> singularity(PolyMesh2D::PoissonSingularity(0.5, 0.0, Eigen::Vector2d(0.0, 0.0)));

    // std::cout << "\n\n" << singularity.value(Eigen::Vector2d(1.0, 0.0)) << "\n";
    // std::cout << "\n\n" << singularity.value(Eigen::Vector2d(0.0, 0.0)) << "\n";
    // std::cout << "\n\n" << singularity.value(Eigen::Vector2d(-1.0, 0.0)) << "\n";
    // std::cout << "\n\n" << singularity.value(Eigen::Vector2d(0.0, 1.0)) << "\n";

    // double osc = Math::PI / 4.0;

    // std::function<double(Eigen::Vector2d)> harmonic_exp_trig = [osc](Eigen::Vector2d x)
    // -> double {return std::exp(osc * x(0)) * std::sin(osc * x(1));};
    // // -> double {return std::exp(osc * x(0)) * std::cos(osc * x(1));};

    // std::function<Eigen::RowVector2d(Eigen::Vector2d)> harmonic_exp_trig_deriv = [osc](Eigen::Vector2d x)
    // -> Eigen::RowVector2d {return osc * std::exp(osc * x(0)) * Eigen::RowVector2d(std::sin(osc * x(1)), std::cos(osc * x(1)));};
    // // -> Eigen::RowVector2d {return osc * std::exp(osc * x(0)) * Eigen::RowVector2d(std::cos(osc * x(1)), -std::sin(osc * x(1)));};

    // PolyMesh2D::Functional::Function<2, 1> singularity(harmonic_exp_trig, harmonic_exp_trig_deriv);

    // std::function<double(Eigen::Vector2d)> harmonic_poly = [](Eigen::Vector2d x) -> double {return std::pow(x(0), 4) - 6.0 * std::pow(x(0), 2) * std::pow(x(1), 2) + std::pow(x(1), 4);};
    // std::function<Eigen::RowVector2d(Eigen::Vector2d)> harmonic_poly_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d {return Eigen::RowVector2d(4.0 * std::pow(x(0), 3) - 12.0 * x(0) * std::pow(x(1), 2), 4.0 * std::pow(x(1), 3) - 12.0 * x(1) * std::pow(x(0), 2));};

    // std::function<double(const Eigen::Vector2d &x)> source_lam = [](const Eigen::Vector2d &x) -> double {return 0.0;};

    // std::function<double(const Eigen::Vector2d &x)> source_lam = [](const Eigen::Vector2d &x) -> double {return -4.0;};

    // std::function<double(const Eigen::Vector2d &x)> poly_zero_bc = [](const Eigen::Vector2d &x) -> double {return x(0) * x(1) * (1.0 - x(0)) * (1.0 - x(1));};
    // std::function<Eigen::RowVector2d(const Eigen::Vector2d &x)> poly_zero_bc_deriv = [](const Eigen::Vector2d &x) -> Eigen::RowVector2d {return Eigen::RowVector2d(x(1) * (1.0 - x(1)) * (1.0 - 2.0 * x(0)), x(0) * (1.0 - x(0)) * (1.0 - 2.0 * x(1)));};

    // std::function<double(const Eigen::Vector2d &x)> source_poly_zero_bc = [](const Eigen::Vector2d &x) -> double {return 2.0 * (x(0) * (1.0 - x(0)) + x(1) * (1.0 - x(1)));};

    std::function<double(Eigen::Vector2d)> reg_val = [](Eigen::Vector2d x) -> double
    { return std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)); };

    std::function<Eigen::RowVector2d(Eigen::Vector2d)> reg_deriv = [](Eigen::Vector2d x) -> Eigen::RowVector2d
    { return Math::PI * Eigen::RowVector2d(std::cos(Math::PI * x(0)) * std::sin(Math::PI * x(1)), std::sin(Math::PI * x(0)) * std::cos(Math::PI * x(1))); };

    std::function<double(const Eigen::Vector2d &x)> source_lam = [](const Eigen::Vector2d &x) -> double { return 2.0 * Math::PI * Math::PI * std::sin(Math::PI * x(0)) * std::sin(Math::PI * x(1)); };

    ScalarFunction2D source(source_lam);
    ScalarFunction2D regularity(reg_val, reg_deriv);
    // ScalarFunction2D regularity(harmonic_exp_trig, harmonic_exp_trig_deriv);
    // ScalarFunction2D exact_sol(poly_zero_bc, poly_zero_bc_deriv);
    ScalarFunction2D enrich(singularity);
    // ScalarFunction2D enrich2(regularity);

    double enrich_tol = 0.15;

    if (false)
    {
        for (size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
        {
            if(locally_enrich(curved_mesh->cell(iT), enrich_tol))
            {
                vem.enrich_cell_basis(iT, enrich);
            }
        }

        for (size_t iE = 0; iE < curved_mesh->n_edges(); ++iE)
        {
            if(locally_enrich(curved_mesh->edge(iE), enrich_tol))
            {
                vem.enrich_edge_basis(iE, PolyMesh2D::Functional::trace(enrich, curved_mesh->edge(iE)->parameterisation()));
            }
        }
    }

    Model model(vem, source);

    model.assemble(params.use_threads);

    Eigen::VectorXd interp_Uvec(vem.interpolate(singularity) + vem.interpolate(regularity));
    // Eigen::VectorXd interp_Uvec(vem.interpolate(singularity));
    // Eigen::VectorXd interp_Uvec(vem.interpolate(exact_sol));

    Eigen::VectorXd UDir = Eigen::VectorXd::Zero(interp_Uvec.size());

    size_t fixed_dofs = 0;

    for (size_t ibF = 0; ibF < curved_mesh->n_b_edges(); ++ibF)
    {
        const size_t iF = curved_mesh->b_edge(ibF)->global_index();
        const size_t offset = vem.global_offset_E(iF);
        for (size_t i = 0; i < vem.local_low_edge_dofs(iF); ++i)
        {
            UDir(offset + i) = interp_Uvec(offset + i);
            fixed_dofs++;
        }
    }
    for (size_t ibV = 0; ibV < curved_mesh->n_b_vertices(); ++ibV)
    {
        const size_t iV = curved_mesh->b_vertex(ibV)->global_index();
        const size_t offset = vem.global_offset_V(iV);

        UDir(offset) = interp_Uvec(offset);
        fixed_dofs++;
    }

    Eigen::VectorXd approx_Uvec(model.solve(UDir));

    std::vector<Eigen::VectorXd> elliptic_projectors;
    elliptic_projectors.reserve(curved_mesh->n_cells());
    for (size_t iT = 0; iT < curved_mesh->n_cells(); ++iT)
    {
        Eigen::VectorXd loc_projector = vem.elliptic_projection(iT, singularity) + vem.elliptic_projection(iT, regularity);
        // Eigen::VectorXd loc_projector = vem.elliptic_projection(iT, singularity);
        elliptic_projectors.push_back(loc_projector);
    }
    auto errors = model.compute_errors(approx_Uvec, interp_Uvec, elliptic_projectors);

    double l2_error = errors[0];
    double h1_error = errors[1];
    double energy_error = errors[2];

    std::cout << "     L2 Error = " << l2_error << "\n";
    std::cout << "     H1 Error = " << h1_error << "\n";
    std::cout << "     Energy Error = " << energy_error << "\n";

    // model.plot(approx_Uvec, interp_Uvec, params.plot_file, singularity + regularity); // elliptic plot and diff plot will not be accurate

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
    out << "DOFs: " << vem.total_edge_dofs() + curved_mesh->n_vertices() - fixed_dofs << "\n";
    out << "NbVertices: " << curved_mesh->n_vertices() << "\n";
    out.close();

    return 0;
}

ModelParameters::ModelParameters(const int argc, const char **argv)
{
    namespace po = boost::program_options;

    // Program options
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Produce help message")("mesh,m", po::value<std::string>(), "Set the mesh")("celldegree,l", po::value<int>(), "Set the degree of the cell polynomials")("edgedegree,k", po::value<int>(), "Set the degree of the edge polynomials")("plot,p", po::value<std::string>(), "Plot to file")("use_threads,u", po::value<bool>(), "Using multithreading")("orthonormalise,o", po::value<bool>(), "Orthonormalise the basis functions");

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
    mesh_name = (vm.count("mesh") ? vm["mesh"].as<std::string>() : "../../../typ2_meshes/mesh2_1.typ2");

    // Get polynomial degrees
    edge_degree = (vm.count("edgedegree") ? vm["edgedegree"].as<int>() : 2);
    cell_degree = std::max(0, edge_degree - 2);

    // Get plot file
    plot_file = (vm.count("plot") ? vm["plot"].as<std::string>() : "");

    // Get use_threads
    use_threads = (vm.count("use_threads") ? vm["use_threads"].as<bool>() : true);

    orthonormalise = (vm.count("orthonormalise") ? vm["orthonormalise"].as<bool>() : true);
}

Model::Model(const XVEMCore &vem, const ScalarFunction2D &src) : m_vem(vem), m_src(src), mesh_ptr(m_vem.get_mesh()) {}

void Model::assemble(bool use_threads)
{
    std::cout << "\n[Scheme] Assembling the linear system\n";

    // Start the timer
    boost::timer::cpu_timer assembly_timer;
    assembly_timer.start();

    // Set up triplets for sparse matrix initialisation
    std::vector<Eigen::Triplet<double>> triplets_GlobMat;

    GlobRHS = Eigen::VectorXd::Zero(m_vem.total_dofs());

    GlobMat.resize(m_vem.total_dofs(), m_vem.total_dofs());

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
            local_poisson_operator(iT, AT[iT], PT[iT]);
            local_source_term(iT, bT[iT]);
        }
    };

    // Running the local constructions in parallel
    parallel_for(mesh_ptr->n_cells(), construct_all_local_contributions, use_threads);

    // Set the local matrices into the correct positions of the global matrices
    for (size_t iT = 0; iT < mesh_ptr->n_cells(); iT++)
    {
        // cell_dofs
        for (size_t i = 0; i < m_vem.local_low_cell_dofs(iT); ++i)
        {
            const size_t iGlobal = m_vem.global_offset_T(iT) + i;
            const size_t iLocal = m_vem.local_offset_T(iT) + i;

            GlobRHS(iGlobal) = bT[iT](iLocal);
            for (size_t j = 0; j < m_vem.local_low_cell_dofs(iT); ++j)
            {
                size_t jGlobal = m_vem.global_offset_T(iT) + j;
                size_t jLocal = m_vem.local_offset_T(iT) + j;
                triplets_GlobMat.emplace_back(iGlobal, jGlobal, AT[iT](iLocal, jLocal));
            }
        }

        PolyMesh2D::CurvedMesh::Cell *cell = mesh_ptr->cell(iT);
        for (size_t iTE = 0; iTE < cell->n_edges(); iTE++)
        {
            const size_t iE = cell->edge(iTE)->global_index();
            for (size_t i = 0; i < m_vem.local_low_edge_dofs(iE); ++i)
            {
                const size_t iGlobal = m_vem.global_offset_E(iE) + i;
                const size_t iLocal = m_vem.local_offset_E(iT, iTE) + i;

                // cell_edge_dofs, edge_cell_dofs
                for (size_t j = 0; j < m_vem.local_low_cell_dofs(iT); ++j)
                {
                    const size_t jGlobal = m_vem.global_offset_T(iT) + j;
                    const size_t jLocal = m_vem.local_offset_T(iT) + j;
                    triplets_GlobMat.emplace_back(iGlobal, jGlobal, AT[iT](iLocal, jLocal));
                    triplets_GlobMat.emplace_back(jGlobal, iGlobal, AT[iT](jLocal, iLocal));
                }

                // edge_edge_dofs
                for (size_t jTE = 0; jTE < cell->n_edges(); jTE++)
                {
                    const size_t jE = cell->edge(jTE)->global_index();
                    for (size_t j = 0; j < m_vem.local_low_edge_dofs(jE); ++j)
                    {
                        const size_t jGlobal = m_vem.global_offset_E(jE) + j;
                        const size_t jLocal = m_vem.local_offset_E(iT, jTE) + j;
                        triplets_GlobMat.emplace_back(iGlobal, jGlobal, AT[iT](iLocal, jLocal));
                    }
                }

                // edge_vertex_dofs, vertex_edge_dofs
                for (size_t jTV = 0; jTV < cell->n_vertices(); jTV++)
                {
                    const size_t jGlobal = m_vem.global_offset_V(cell->vertex(jTV)->global_index());
                    const size_t jLocal = m_vem.local_offset_V(iT, jTV);
                    triplets_GlobMat.emplace_back(iGlobal, jGlobal, AT[iT](iLocal, jLocal));
                    triplets_GlobMat.emplace_back(jGlobal, iGlobal, AT[iT](jLocal, iLocal));
                }
            }
        }

        for (size_t iTV = 0; iTV < cell->n_vertices(); iTV++)
        {
            const size_t iGlobal = m_vem.global_offset_V(cell->vertex(iTV)->global_index());
            const size_t iLocal = m_vem.local_offset_V(iT, iTV);

            // cell_vertex_dofs, vertex_cell_dofs
            for (size_t j = 0; j < m_vem.local_low_cell_dofs(iT); ++j)
            {
                const size_t jGlobal = m_vem.global_offset_T(iT) + j;
                const size_t jLocal = m_vem.local_offset_T(iT) + j;
                triplets_GlobMat.emplace_back(iGlobal, jGlobal, AT[iT](iLocal, jLocal));
                triplets_GlobMat.emplace_back(jGlobal, iGlobal, AT[iT](jLocal, iLocal));
            }

            // vertex_vertex_dofs
            for (size_t jTV = 0; jTV < cell->n_vertices(); jTV++)
            {
                const size_t jGlobal = m_vem.global_offset_V(cell->vertex(jTV)->global_index());
                const size_t jLocal = m_vem.local_offset_V(iT, jTV);
                triplets_GlobMat.emplace_back(iGlobal, jGlobal, AT[iT](iLocal, jLocal));
            }
        }
    }

    // Construct the global matrices
    GlobMat.setFromTriplets(std::begin(triplets_GlobMat), std::end(triplets_GlobMat));

    std::cout << "     Assembly time = " << assembly_timer.elapsed().wall * pow(10, -9) << "s\n";
}

Eigen::VectorXd Model::solve(const Eigen::VectorXd &UDir)
{
    std::cout << "\n[Scheme] Solving the linear system\n";

    // Start the timer
    boost::timer::cpu_timer solve_timer;
    solve_timer.start();

    Eigen::SparseMatrix<double> SysMat = GlobMat;
    for (size_t ibF = 0; ibF < mesh_ptr->n_b_edges(); ++ibF)
    {
        const size_t iF = mesh_ptr->b_edge(ibF)->global_index();
        const int n_local_low_edge_dofs = int(m_vem.local_low_edge_dofs(iF));
        const int offset = int(m_vem.global_offset_E(iF));
        for (int i = 0; i < n_local_low_edge_dofs; ++i)
        {
            SysMat.row(offset + i) *= 0.0;
        }
    }
    for (size_t ibV = 0; ibV < mesh_ptr->n_b_vertices(); ++ibV)
    {
        const size_t iV = mesh_ptr->b_vertex(ibV)->global_index();
        const int offset = int(m_vem.global_offset_V(iV));

        SysMat.row(offset) *= 0.0;
    }

    Eigen::VectorXd dirichlet_terms = SysMat * UDir;

    for (size_t ibF = 0; ibF < mesh_ptr->n_b_edges(); ++ibF)
    {
        const size_t iF = mesh_ptr->b_edge(ibF)->global_index();
        const int n_local_low_edge_dofs = int(m_vem.local_low_edge_dofs(iF));
        const int offset = int(m_vem.global_offset_E(iF));
        for (int i = 0; i < n_local_low_edge_dofs; ++i)
        {
            SysMat.col(offset + i) *= 0.0;
            SysMat.coeffRef(offset + i, offset + i) = 1.0;
        }
    }
    for (size_t ibV = 0; ibV < mesh_ptr->n_b_vertices(); ++ibV)
    {
        const size_t iV = mesh_ptr->b_vertex(ibV)->global_index();
        const int offset = int(m_vem.global_offset_V(iV));

        SysMat.col(offset) *= 0.0;
        SysMat.coeffRef(offset, offset) = 1.0;
    }
    SysMat.prune(0.0);

    Eigen::VectorXd B = GlobRHS - dirichlet_terms;

    // Solve the statically condesed system using BiCGSTAB
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(SysMat);

    // xF = INV( AII - AIT * INV_ATT * ATI ) * ( BI - AIT * INV_ATT * BT - (AID - AIT * INV_ATT * ATD * UD) )
    Eigen::VectorXd Xh = solver.solve(B);

    // Print solver iterations and estimated error
    std::cout << "     [solver] #iterations: " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;

    std::cout << "     Solving time = " << solve_timer.elapsed().wall * pow(10, -9) << "s\n";

    return Xh + UDir;
}

std::vector<double> Model::compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::vector<Eigen::VectorXd> &elliptic_projectors)
{
    std::cout << "\n[Scheme] Computing error";
    boost::timer::cpu_timer error_timer;
    error_timer.start();

    // for (size_t iT = 0; iT < mesh_ptr->n_cells(); ++iT)
    // {
    //     Eigen::VectorXd local_approx = m_vem.restr(approx_Uvec, iT);
    //     Eigen::VectorXd local_interp = m_vem.restr(interp_Uvec, iT);

    //     Eigen::VectorXd diff_vec = local_approx - local_interp;

    //     energy_diff += diff_vec.transpose() * AT[iT] * diff_vec;
    //     energy_exact += local_interp.transpose() * AT[iT] * local_interp;

    //     auto quad_handle(m_vem.get_quad_handle());

    //     Eigen::MatrixXd ST = quad_handle.h1_product(*m_vem.high_cell_basis(iT), *m_vem.high_cell_basis(iT), mesh_ptr->cell(iT), true);
    //     Eigen::MatrixXd MT = quad_handle.l2_product(*m_vem.high_cell_basis(iT), *m_vem.high_cell_basis(iT), mesh_ptr->cell(iT), true);

    //     Eigen::VectorXd SigmaT = elliptic_projectors[iT];

    //     Eigen::VectorXd pT_minus_SigmaT = SigmaT - PT[iT] * local_approx;

    //     H1_diff += pT_minus_SigmaT.transpose() * ST * pT_minus_SigmaT;
    //     H1_exact += SigmaT.transpose() * ST * SigmaT;

    //     L2_diff += pT_minus_SigmaT.transpose() * MT * pT_minus_SigmaT;
    //     L2_exact += SigmaT.transpose() * MT * SigmaT;
    // }

    std::vector<double> energy_diff_vec(mesh_ptr->n_cells());
    std::vector<double> energy_exact_vec(mesh_ptr->n_cells());
    std::vector<double> L2_diff_vec(mesh_ptr->n_cells());
    std::vector<double> L2_exact_vec(mesh_ptr->n_cells());
    std::vector<double> H1_diff_vec(mesh_ptr->n_cells());
    std::vector<double> H1_exact_vec(mesh_ptr->n_cells());

    // Create extended bases
    std::function<void(size_t, size_t)> compute_errors_in_parallel = [&](size_t start, size_t end) -> void
    {
        for (size_t iT = start; iT < end; iT++)
        {
            Eigen::VectorXd local_approx = m_vem.restr(approx_Uvec, iT);
            Eigen::VectorXd local_interp = m_vem.restr(interp_Uvec, iT);

            Eigen::VectorXd diff_vec = local_approx - local_interp;

            energy_diff_vec[iT] = diff_vec.transpose() * AT[iT] * diff_vec;
            energy_exact_vec[iT] = local_interp.transpose() * AT[iT] * local_interp;

            auto quad_handle(m_vem.get_quad_handle());

            Eigen::MatrixXd ST = quad_handle.h1_product(*m_vem.high_cell_basis(iT), *m_vem.high_cell_basis(iT), mesh_ptr->cell(iT), true);
            Eigen::MatrixXd MT = quad_handle.l2_product(*m_vem.high_cell_basis(iT), *m_vem.high_cell_basis(iT), mesh_ptr->cell(iT), true);

            Eigen::VectorXd SigmaT = elliptic_projectors[iT];
            Eigen::VectorXd pT_minus_SigmaT = SigmaT - PT[iT] * local_approx;

            H1_diff_vec[iT] = pT_minus_SigmaT.transpose() * ST * pT_minus_SigmaT;
            H1_exact_vec[iT] = SigmaT.transpose() * ST * SigmaT;

            L2_diff_vec[iT] = pT_minus_SigmaT.transpose() * MT * pT_minus_SigmaT;
            L2_exact_vec[iT] = SigmaT.transpose() * MT * SigmaT;
        }
    };
    parallel_for(mesh_ptr->n_cells(), compute_errors_in_parallel, true);

    double energy_diff = std::accumulate(energy_diff_vec.begin(), energy_diff_vec.end(), 0.0);
    double energy_exact = std::accumulate(energy_exact_vec.begin(), energy_exact_vec.end(), 0.0);

    double L2_diff = std::accumulate(L2_diff_vec.begin(), L2_diff_vec.end(), 0.0);
    double L2_exact = std::accumulate(L2_exact_vec.begin(), L2_exact_vec.end(), 0.0);

    double H1_diff = std::accumulate(H1_diff_vec.begin(), H1_diff_vec.end(), 0.0);
    double H1_exact = std::accumulate(H1_exact_vec.begin(), H1_exact_vec.end(), 0.0);

    error_timer.stop();
    std::cout << " (" << error_timer.elapsed().wall * 1E-9 << "s)\n";

    return {std::sqrt(L2_diff / L2_exact), std::sqrt(H1_diff / H1_exact), std::sqrt(energy_diff / energy_exact)};
}

void Model::plot(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::string &plot_file, const ScalarFunction2D &sol)
{
    if (plot_file != "")
    {
        std::vector<double> exact, potential, elliptic, diff;

        exact.reserve(mesh_ptr->n_vertices());
        potential.reserve(mesh_ptr->n_vertices());
        elliptic.reserve(mesh_ptr->n_vertices());
        diff.reserve(mesh_ptr->n_vertices());

        for (size_t iV = 0; iV < mesh_ptr->n_vertices(); iV++)
        {
            exact.push_back(sol.value(mesh_ptr->vertex(iV)->coords()));

            auto xV = mesh_ptr->vertex(iV)->coords();
            auto cList = mesh_ptr->vertex(iV)->get_cells();

            double loc_potential = 0.0;
            double loc_elliptic = 0.0;
            double loc_diff = 0.0;
            for (size_t ilT = 0; ilT < cList.size(); ilT++)
            {
                size_t iT = cList[ilT]->global_index();

                Eigen::VectorXd local_potential = PT[iT] * m_vem.restr(approx_Uvec, iT);
                Eigen::VectorXd local_elliptic = m_vem.elliptic_projection(iT, sol);
                Eigen::VectorXd local_diff = local_potential - local_elliptic;

                for (size_t i = 0; i < m_vem.local_high_cell_dofs(iT); i++)
                {
                    loc_potential += local_potential(i) * m_vem.high_cell_basis(iT)->value(i, xV);
                    loc_elliptic += local_elliptic(i) * m_vem.high_cell_basis(iT)->value(i, xV);
                    loc_diff += local_diff(i) * m_vem.high_cell_basis(iT)->value(i, xV);
                }
            }
            loc_potential = loc_potential / (cList.size());
            loc_elliptic = loc_elliptic / (cList.size());
            loc_diff = loc_diff / (cList.size());

            // loc_diff = loc_elliptic - exact[iV];

            potential.push_back(loc_potential);
            elliptic.push_back(loc_elliptic);
            diff.push_back(loc_diff);
        }

        VtuWriter plotdata(mesh_ptr);

        plotdata.write_to_vtu("exact-" + plot_file + ".vtu", exact);
        plotdata.write_to_vtu("potential-" + plot_file + ".vtu", potential);
        plotdata.write_to_vtu("elliptic-" + plot_file + ".vtu", elliptic);
        plotdata.write_to_vtu("diff-" + plot_file + ".vtu", diff);
    }
}

void Model::local_poisson_operator(const size_t iT, Eigen::MatrixXd &AT, Eigen::MatrixXd &PT)
{
    CurvedMesh::Cell *cell = mesh_ptr->cell(iT);

    const size_t n_local_low_cell_dofs = m_vem.local_low_cell_dofs(iT);

    auto quad_handle(m_vem.get_quad_handle());

    // Calculate stiffness matrix
    Eigen::MatrixXd ST(quad_handle.h1_product(*m_vem.high_cell_basis(iT), *m_vem.high_cell_basis(iT), cell, true));
    Eigen::MatrixXd M_LOWT_LOWT = quad_handle.l2_product(*m_vem.low_cell_basis(iT), *m_vem.low_cell_basis(iT), cell, true);
    Eigen::MatrixXd M_LOWT_HIGHT = quad_handle.l2_product(*m_vem.low_cell_basis(iT), *m_vem.high_cell_basis(iT), cell, false);

    PT = m_vem.potential_reconstruction(iT);
    Eigen::MatrixXd ATCONS = PT.transpose() * ST * PT;

    Eigen::MatrixXd DT = M_LOWT_LOWT.inverse() * M_LOWT_HIGHT * PT;
    DT.topLeftCorner(n_local_low_cell_dofs, n_local_low_cell_dofs) -= Eigen::MatrixXd::Identity(n_local_low_cell_dofs, n_local_low_cell_dofs);

    double hT = cell->diam();

    Eigen::MatrixXd ATSTAB = std::pow(hT, -2) * DT.transpose() * M_LOWT_LOWT * DT;

    for (size_t iTE = 0; iTE < cell->n_edges(); iTE++)
    {
        CurvedMesh::Edge *edge = cell->edge(iTE);

        const size_t iE = edge->global_index();

        Eigen::MatrixXd M_HIGHE_HIGHE = quad_handle.l2_product(*m_vem.high_edge_basis(iE), *m_vem.high_edge_basis(iE), edge, true);

        Eigen::MatrixXd M_HIGHE_HIGHT = quad_handle.l2_product(*m_vem.high_edge_basis(iE), *m_vem.high_cell_basis(iT), edge, false);

        Eigen::MatrixXd DTE = M_HIGHE_HIGHE.inverse() * M_HIGHE_HIGHT * PT - m_vem.continuous_reconstruction(iT, iTE);

        ATSTAB += std::pow(hT, -1) * DTE.transpose() * M_HIGHE_HIGHE * DTE;
    }

    AT = ATCONS + ATSTAB;
}

void Model::local_source_term(const size_t iT, Eigen::VectorXd &bT)
{
    bT = Eigen::VectorXd::Zero(m_vem.local_dofs(iT));

    auto quad_handle(m_vem.get_quad_handle());
    bT.head(m_vem.local_low_cell_dofs(iT)) = quad_handle.l2_product(m_src, *m_vem.low_cell_basis(iT), mesh_ptr->cell(iT));
}
