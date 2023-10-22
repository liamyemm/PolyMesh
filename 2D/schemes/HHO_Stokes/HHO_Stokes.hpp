#include "Mesh.hpp"
#include "StokesCore.hpp"
#include "StokesTests.hpp"

#include <unsupported/Eigen/SparseExtra> // Eigen::SparseMatrix
#include <Eigen/Dense>                   // Eigen::Matrix

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#endif

#ifndef HHO_STOKES_HPP
#define HHO_STOKES_HPP

namespace PolyMesh2D
{
    namespace HHOSTOKES
    {
        // using Quadrature::GaussLegendre1D;
        using Quadrature::QuadratureRule;

        using CurvedMesh::Cell;
        using CurvedMesh::Edge;
        using CurvedMesh::Mesh;
        using CurvedMesh::Vertex;

        using Functional::ScalarFunction2D;
        using Functional::VectorFunction2D;

        using Functional::Curve;

        struct ModelParameters
        {
            ModelParameters(const int argc, const char **argv); ///< constructor

            std::string mesh_name, plot_file;
            bool use_threads, orthonormalise;
            unsigned cell_degree, edge_degree;
            double radius;
        };
        class Model
        {
        public:
            Model(const StokesCore &stokes_backend, const VectorFunction2D &src);
            void assemble(bool use_threads = true);
            Eigen::VectorXd solve(const Eigen::VectorXd &UDir);
            std::vector<double> compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::vector<Eigen::VectorXd> &elliptic_projectors, const std::vector<Eigen::VectorXd> &l2_velocity_projectors);
            void plot(const Eigen::VectorXd &approx_Uvec);

        private:
            void local_stokes_operator(const size_t iT, Eigen::MatrixXd &AT, Eigen::MatrixXd &RT);
            void local_source_term(const size_t iT, Eigen::VectorXd &bT);

            const StokesCore &m_stokes;
            const VectorFunction2D &m_src;

            Mesh *mesh_ptr;

            std::vector<Eigen::MatrixXd> RT;
            std::vector<Eigen::MatrixXd> DivT;
            std::vector<Eigen::MatrixXd> AT;

            Eigen::SparseMatrix<double> inv_LTll_LTlg;
            Eigen::SparseMatrix<double> LTgl;
            Eigen::SparseMatrix<double> LTgg;

            Eigen::VectorXd inv_LTll_rTl;
            Eigen::VectorXd rTg;
        };
    }
}

#endif