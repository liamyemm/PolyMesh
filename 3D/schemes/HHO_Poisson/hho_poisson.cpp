#include "hho_poisson.hpp"
#include "mesh_builder.hpp"

using namespace PolyMesh3D;

int main(const int argc, const char **argv)
{
    // Build mesh and reorder edges
    StraightMesh::MeshBuilder builder(argv[1]);
    std::unique_ptr<StraightMesh::Mesh> mesh;

    try
    {
        mesh = builder.build_the_mesh();
    }
    catch (std::string msg)
    {
        std::cerr << msg;
    }

    bool valid = mesh->test();

    std::cout << "\n";

    Quadrature::GaussLegendre1D quad(10);
    // Quadrature::GaussLegendre1D high_quad(30);

    std::vector<Quadrature::QuadratureRule> edge_quads;
    std::vector<Quadrature::QuadratureRule> face_quads;
    std::vector<Quadrature::QuadratureRule> cell_quads;

    for (auto &E : mesh->get_edges())
    {
        edge_quads.push_back(generate_quadrature_rule(E, quad));
    }

    for (auto &F : mesh->get_faces())
    {
        // std::vector<QuadratureRule<double>> local_edge_quads(c->n_edges());
        std::vector<Quadrature::QuadratureRule> local_edge_quads;
        for (size_t iFE = 0; iFE < F->n_edges(); ++iFE)
        {
            local_edge_quads.push_back(edge_quads[F->edge(iFE)->global_index()]);
        }
        face_quads.push_back(generate_quadrature_rule(F, local_edge_quads, quad));
    }

    for (auto &C : mesh->get_cells())
    {
        // std::vector<QuadratureRule<double>> local_edge_quads(c->n_edges());
        std::vector<Quadrature::QuadratureRule> local_face_quads;
        for (size_t iTF = 0; iTF < C->n_faces(); ++iTF)
        {
            local_face_quads.push_back(face_quads[C->face(iTF)->global_index()]);
        }
        cell_quads.push_back(generate_quadrature_rule(C, local_face_quads, quad));
    }



    for (size_t iE = 0; iE < mesh->n_edges(); ++iE)
    {
        double measure_integral = 0.0;
        StraightMesh::VectorRd center_integral = StraightMesh::VectorRd::Zero();
        for (auto & quad_i : edge_quads[iE])
        {
            measure_integral += quad_i.w;
            center_integral += quad_i.w * quad_i.x;
        }
        center_integral /= measure_integral;

        if(std::abs(measure_integral - mesh->edge(iE)->measure()) > 1E-12)
        {
            std::cerr << "Discrepancy in measures: edge measure is " << mesh->edge(iE)->measure() << " integral gives " << measure_integral << "\n";
            valid = false;
        }

        if((center_integral - mesh->edge(iE)->center_mass()).norm() > 1E-12)
        {
            std::cerr << "Discrepancy in center masses: edge center is (" << mesh->edge(iE)->center_mass().transpose() << ") integral gives (" << center_integral.transpose() << ")\n";
            valid = false;
        }

    }

    for (size_t iF = 0; iF < mesh->n_faces(); ++iF)
    {
        double measure_integral = 0.0;
        StraightMesh::VectorRd center_integral = StraightMesh::VectorRd::Zero();
        for (auto & quad_i : face_quads[iF])
        {
            measure_integral += quad_i.w;
            center_integral += quad_i.w * quad_i.x;
        }
        center_integral /= measure_integral;

        if(std::abs(measure_integral - mesh->face(iF)->measure()) > 1E-12)
        {
            std::cerr << "Discrepancy in measures: face measure is " << mesh->face(iF)->measure() << " integral gives " << measure_integral << "\n";
            valid = false;
        }

        if((center_integral - mesh->face(iF)->center_mass()).norm() > 1E-12)
        {
            std::cerr << "Discrepancy in center masses: face center is (" << mesh->face(iF)->center_mass().transpose() << ") integral gives (" << center_integral.transpose() << ")\n";
            valid = false;
        }

    }

    double total_volume = 0.0;
    StraightMesh::VectorRd total_center = StraightMesh::VectorRd::Zero();

    for (size_t iT = 0; iT < mesh->n_cells(); ++iT)
    {
        double measure_integral = 0.0;
        StraightMesh::VectorRd center_integral = StraightMesh::VectorRd::Zero();
        for (auto & quad_i : cell_quads[iT])
        {
            measure_integral += quad_i.w;
            center_integral += quad_i.w * quad_i.x;
            total_center += quad_i.w * quad_i.x;
        }
        center_integral /= measure_integral;

        total_volume += measure_integral;

        if(std::abs(measure_integral - mesh->cell(iT)->measure()) > 1E-12)
        {
            std::cerr << "Discrepancy in measures: cell measure is " << mesh->cell(iT)->measure() << " integral gives " << measure_integral << "\n";
            valid = false;
        }

        if((center_integral - mesh->cell(iT)->center_mass()).norm() > 1E-12)
        {
            std::cerr << "Discrepancy in center masses: cell center is (" << mesh->cell(iT)->center_mass().transpose() << ") integral gives (" << center_integral.transpose() << ")\n";
            valid = false;
        }

    }

    total_center /= total_volume;

    if(std::abs(total_volume - 1.0) > 1E-12)
    {
        std::cerr << "Discrepancy in total volume: integral gives " << total_volume << "\n";
        valid = false;
    }

    if((total_center - StraightMesh::VectorRd(0.5, 0.5, 0.5)).norm() > 1E-12)
    {
        std::cerr << "Discrepancy in total center: integral gives (" << total_center.transpose() << ")\n";
        valid = false;
    }

    return !valid;
}