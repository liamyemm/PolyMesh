#include "Mesh.hpp"
#include "XVEMCore.h"

#include <unsupported/Eigen/SparseExtra> // Eigen::SparseMatrix
#include <Eigen/Dense> // Eigen::Matrix

namespace PolyMesh2D
{
    namespace XVEM
    {
        using CurvedMesh::Mesh;

        using Functional::ScalarFunction2D;

        struct ModelParameters
        {
            ModelParameters(const int argc, const char **argv); ///< constructor

            std::string mesh_name, plot_file;
            bool use_threads, orthonormalise;
            int cell_degree, edge_degree;
        };

        class Model
        {
        public:
            Model(const XVEMCore &vem, const ScalarFunction2D &src);
            void assemble(bool use_threads = true);
            Eigen::VectorXd solve(const Eigen::VectorXd &UDir);
            std::vector<double> compute_errors(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::vector<Eigen::VectorXd> &elliptic_projectors);
            void plot(const Eigen::VectorXd &approx_Uvec, const Eigen::VectorXd &interp_Uvec, const std::string &plot_file, const ScalarFunction2D &sol);

        private:
            void local_poisson_operator(const size_t iT, Eigen::MatrixXd &AT, Eigen::MatrixXd &PT);
            void local_source_term(const size_t iT, Eigen::VectorXd &bT);

            const XVEMCore &m_vem;
            const ScalarFunction2D &m_src;

            Mesh *mesh_ptr;

            std::vector<Eigen::MatrixXd> PT;
            std::vector<Eigen::MatrixXd> AT;

            Eigen::VectorXd GlobRHS;
            Eigen::SparseMatrix<double> GlobMat;
        };
    }
}