// Class to provide various test cases for the homogeneous Poisson problem on the unit circle
//
//
// Author: Liam Yemm (liam.yemm@monash.edu)
//

#ifndef _TEST_CASE_HPP
#define _TEST_CASE_HPP

#include "function.hpp"
#include "Mesh.hpp"
#include "math.hpp"

/*!
 * @defgroup TestCases
 *	@brief Defines test cases (exact solutions, source terms)
 */

namespace PolyMesh2D
{
    namespace HHOPOISSON
    {
        using PolyMesh2D::CurvedMesh::Mesh;
        using Functional::ScalarFunction2D;
        using Functional::Curve;

        // @addtogroup TestCases
        //@{

        // ----------------------------------------------------------------------------
        //                            Class definition
        // ----------------------------------------------------------------------------

        /// The TestCase class provides definition of test cases
        class TestCase
        {
        public:
            typedef typename ScalarFunction2D::InputType InputType;
            typedef typename ScalarFunction2D::OutputType OutputType;
            typedef typename ScalarFunction2D::DerivativeType DerivativeType;

            /// Initialise data
            TestCase(
                unsigned id_test_case, ///< The id of the test case.
                char bdry = 'C'///< The id of the bdry curve: 'C' = circle, 'E' = ellipse
            );

            /// Returns the exact solution
            ScalarFunction2D sol();

            /// Returns the exact source term
            ScalarFunction2D src();

            Curve get_boundary_param();
            ScalarFunction2D get_level_set();

        private:
            unsigned m_id_test_case;

            Curve bdry_param;
            ScalarFunction2D level_set;

            std::function<double(Functional::ColVector)> level_set_laplace;
            std::function<double(double)> u = [](const double t) -> double
            {
                return 0.0;
            };
            std::function<double(double)> Du = u;
            std::function<double(double)> DDu = u;
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

#endif //_TEST_CASE_HPP
