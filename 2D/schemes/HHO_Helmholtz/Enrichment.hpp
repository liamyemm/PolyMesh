// Class to provide various test cases for the homogeneous Poisson problem on the unit circle
//
//
// Author: Liam Yemm (liam.yemm@monash.edu)
//

#ifndef _ENRICHMENT_HPP
#define _ENRICHMENT_HPP

#include "function.hpp"
#include "HybridCore.hpp"
#include "Mesh.hpp"
#include "math.hpp"

namespace PolyMesh2D
{
    namespace HHOHELMHOLTZ
    {
        using Functional::Curve;
        using Functional::ScalarFunction1D;
        using Functional::ScalarFunction2D;
        using PolyMesh2D::CurvedMesh::Mesh;

        // @addtogroup TestCases
        //@{

        // ----------------------------------------------------------------------------
        //                            Class definition
        // ----------------------------------------------------------------------------

        /// The TestCase class provides definition of test cases
        class Enrichment
        {
        public:
            typedef typename ScalarFunction2D::InputType InputType;
            typedef typename ScalarFunction2D::OutputType OutputType;
            typedef typename ScalarFunction2D::DerivativeType DerivativeType;

            /// Initialise data
            Enrichment(
                double lambda);

            void globally_enrich(HybridCore &hho);
            void locally_enrich(HybridCore &hho);

            // ScalarFunction2D loc_enrichment(const Eigen::Vector2d &cell_center) const;
            // ScalarFunction2D loc_enrichment_laplace(const Eigen::Vector2d &cell_center) const;
            // ScalarFunction1D loc_enrichment_neumann_trace(const Eigen::Vector2d &cell_center, const Curve &edge_param) const;

        private:
            const double m_lambda;

            ScalarFunction2D glob_enrichment() const;
            ScalarFunction2D glob_enrichment_laplace() const;
            ScalarFunction1D glob_enrichment_neumann_trace(const Curve &edge_param) const;

            ScalarFunction2D loc_enrichment(const Eigen::Vector2d &cell_center) const;
            ScalarFunction2D loc_enrichment_laplace(const Eigen::Vector2d &cell_center) const;
            ScalarFunction1D loc_enrichment_neumann_trace(const Eigen::Vector2d &cell_center, const Curve &edge_param) const;
        };

        inline void reorder_edges(Mesh *mesh_ptr)
        {
            // Create vector with all non-Dirichlet edges first, and all Dirichlet edges at the end
            std::vector<size_t> new_to_old(mesh_ptr->n_edges(), 0);
            // Index for non-Dirichlet edges start at 0 and increases, index for Dirichlet edges start at n_edges()-1
            // and decreases
            size_t idx_nondir = 0;
            size_t idx_dir = mesh_ptr->n_edges() - 1;
            for (auto &edge : mesh_ptr->get_edges())
            {
                assert(idx_nondir < idx_dir + 1);
                if (edge->is_boundary())
                {
                    new_to_old[idx_dir] = edge->global_index();
                    idx_dir--;
                }
                else
                {
                    new_to_old[idx_nondir] = edge->global_index();
                    idx_nondir++;
                }
            }

            // Check: idx_dir and idx_nondir must just have crossed
            assert(idx_nondir == idx_dir + 1);

            // Reordering
            mesh_ptr->renum('E', new_to_old);
        }

        //@}
    }
}

#endif //_ENRICHMENT_HPP
