#include "Mesh.hpp"
#include "HybridCore.hpp"

#include <unsupported/Eigen/SparseExtra> // Eigen::SparseMatrix
#include <Eigen/Dense> // Eigen::Matrix

#ifdef WITH_MKL
#include <Eigen/PardisoSupport>
#endif

namespace PolyMesh2D
{
    namespace HHOHELMHOLTZ
    {
        using Quadrature::GaussLegendre1D;
        using Quadrature::QuadratureRule;

        using CurvedMesh::Mesh;
        using CurvedMesh::Cell;
        using CurvedMesh::Edge;
        using CurvedMesh::Vertex;

        using Functional::ScalarFunction2D;

        using Functional::Curve;

        struct ModelParameters
        {
            ModelParameters(const int argc, const char **argv); ///< constructor

            std::string mesh_name, plot_file;
            bool use_threads, orthonormalise;
            unsigned cell_degree, edge_degree;
            double wave_number;
        };

        class Model
        {
        public:
            Model(const HybridCore &hho, const ScalarFunction2D &src, const double wave_number);
            void assemble(bool use_threads = true);
            Eigen::VectorXd solve(const Eigen::VectorXd &UDir);
            std::vector<double> compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const ScalarFunction2D &sol);
            void plot(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::string &plot_file, const ScalarFunction2D &sol);

        private:
            void local_helmholtz_operator(const size_t iT, Eigen::MatrixXd &AT, Eigen::MatrixXd &PT);
            void local_source_term(const size_t iT, Eigen::VectorXd &bT);

            const HybridCore &m_hho;
            const ScalarFunction2D &m_src;

            const double m_wave_number;

            Mesh *mesh_ptr;

            std::vector<Eigen::MatrixXd> PT;
            std::vector<Eigen::MatrixXd> AT;

            Eigen::VectorXd GlobRHS;
            Eigen::VectorXd ScRHS;

            Eigen::SparseMatrix<double> GlobMat;
            Eigen::SparseMatrix<double> ScBeMat;
        };
    }
}