// Class to read a typ2 mesh file
//
// Author: Liam Yemm (liam.yemm@monash.edu)
//

#ifndef MESHREADERTYP2_HPP
#define MESHREADERTYP2_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>

namespace PolyMesh2D
{
    namespace StraightMesh
    {

        // ----------------------------------------------------------------------------
        //                            Class definition
        // ----------------------------------------------------------------------------

        /*!
         * \addtogroup Mesh
         * @{
         */

        /// The MeshReaderTyp2 class provides functions to read a typ2 mesh file
        class MeshReaderTyp2
        {
        public:
            /**
             * Constructor for mesh reader
             *
             * @param file_name name of the file name, needs to include the full path
             **/
            MeshReaderTyp2(std::string file_name);

            /**
             * Reads the file into the specified containers
             *
             * @param vertices reference to a vector to hold the vertices coordinates
             * @param cells reference to a vector to hold the cell vertices
             **/
            void read_mesh(std::vector<std::array<double, 2>> &vertices, std::vector<std::vector<std::size_t>> &cells);

        private:
            std::string _file_name; ///< name of the file being read
        };

        //@}

    }
} 
#endif
