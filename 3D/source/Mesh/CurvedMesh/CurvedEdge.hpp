#include "CurvedPolytope.hpp"
#include "function.hpp"

#ifndef _CURVEDEDGE_HPP
#define _CURVEDEDGE_HPP

namespace PolyMesh2D
{
    namespace CurvedMesh
    {
        /*!
         * \addtogroup Mesh
         * @{
         */

        class Edge : public Polytope
        {
        public:
            Edge(size_t index, Functional::Curve &parameterisation);

            int vertex_orientation(const size_t vertex_index) const; ///< Return the orientation of the Vertex located at vertex_index

            VectorRd tangent(double t) const;                  ///< Return the tangent of the Edge at the parameter @arg t
            VectorRd normal(double t) const;                   ///< Return the normal of the Edge at the parameter @arg t
            const Functional::Curve &parameterisation() const; ///< Return a const reference to the parameterisation of the Edge

            std::array<VectorRd, 2> coords() const;

            void plot(std::ofstream *out, int partitions) const;

            bool is_straight() const;
            void set_straight();

            bool test() const;

        private:
            Functional::Curve _parameterisation;
            bool _is_straight;
        };

        //@}
    }
}

#endif