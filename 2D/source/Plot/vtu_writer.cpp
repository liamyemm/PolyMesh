// Creates a .vtu file of the solution, to be visualised by paraview
//
// Author: Jerome Droniou (jerome.droniou@monash.edu)
//

// #include "mesh.hpp"
// #include "vertex.hpp"
#include <stdio.h>
#include <vtu_writer.hpp>

// using namespace HArDCore2D;

VtuWriter::VtuWriter(const PolyMesh2D::CurvedMesh::Mesh *mesh) : _mesh(mesh) {}

void VtuWriter::write_header(FILE *pFile)
{
    fprintf(pFile, "%s", "<VTKFile type=\"UnstructuredGrid\"");
    fprintf(pFile, "%s", " version=\"0.1\"");
    fprintf(pFile, "%s\n", " byte_order=\"LittleEndian\">");
    fprintf(pFile, "%s\n", "   <UnstructuredGrid>");
    fprintf(pFile, "%s", "      <Piece  ");
    fprintf(pFile, "%s%i%s", "NumberOfPoints=\"", (int)_mesh->n_vertices(), "\"  ");
    fprintf(pFile, "%s%i%s\n", "NumberOfCells=\"", (int)_mesh->n_cells(), "\">");
}

void VtuWriter::write_vertices(FILE *pFile)
{
    // VERTICES
    fprintf(pFile, "%s\n", "         <Points>");
    fprintf(pFile, "%s%s%s\n",
            "            <DataArray type=\"Float32\" NumberOfComponents=\"",
            "3", "\" format=\"ascii\">");
    for (auto &v : _mesh->get_vertices())
    {
        // vertices coords (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1}
        // z_{i=1}
        //....
        // x_{i=n} y_{i=n} z_{i=n}) i is the node index
        Eigen::Vector2d v_coords = v->coords();
        fprintf(pFile, "%e %e %e", v_coords(0), v_coords(1), 0.0);
        fprintf(pFile, "%s\n", "");
    }

    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "         </Points>");
}

void VtuWriter::write_vertices(FILE *pFile, const std::vector<double> &data)
{
    // VERTICES
    fprintf(pFile, "%s\n", "         <Points>");
    fprintf(pFile, "%s%s%s\n",
            "            <DataArray type=\"Float32\" NumberOfComponents=\"",
            "3", "\" format=\"ascii\">");
    for (size_t iV = 0; iV < _mesh->n_vertices(); iV++)
    {
        // vertices coords (x_{i=0} y_{i=0} z_{i=0} x_{i=1} y_{i=1}
        // z_{i=1}
        //....
        // x_{i=n} y_{i=n} z_{i=n}) i is the node index
        auto v = _mesh->vertex(iV)->coords();
        fprintf(pFile, "%e %e %e", v(0), v(1), data[iV]);
        fprintf(pFile, "%s\n", "");
    }

    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "         </Points>");
}

void VtuWriter::write_cells(FILE *pFile)
{
    fprintf(pFile, "%s\n", "         <Cells>");
    fprintf(pFile, "%s\n",
            "            <DataArray type=\"Int32\" Name=\"connectivity\" "
            "format=\"ascii\">");
    std::vector<PolyMesh2D::CurvedMesh::Cell *> clist = _mesh->get_cells();
    // Lists all vertices of all nodes
    for (auto &cell : clist)
    {
        std::vector<PolyMesh2D::CurvedMesh::Vertex *> vertices = cell->get_vertices();
        for (auto &vertex : vertices)
        {
            fprintf(pFile, "%i ", (int)vertex->global_index());
        }
        fprintf(pFile, "%s\n", "");
    }
    fprintf(pFile, "%s\n", "");
    fprintf(pFile, "%s\n", "            </DataArray>");

    fprintf(pFile, "%s\n",
            "            <DataArray type=\"Int32\" Name=\"offsets\" "
            "format=\"ascii\">");
    size_t prev = 0;
    // probably could be done more efficiently like building fofsets
    // array when making cells list but...
    for (auto &cell : clist)
    {
        prev += cell->n_vertices(); // offsets are the end position of the cell in the nodes array?
        fprintf(pFile, "%i ", (int)prev);
    }
    fprintf(pFile, "%s\n", "");

    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n",
            "            <DataArray type=\"Int32\" Name=\"types\" "
            "format=\"ascii\">");
    // probably could be done more efficiently like building fofsets
    // array when making cells list but...
    for (size_t iC = 0; iC < clist.size(); iC++)
    {
        fprintf(pFile, "%i ", 7);
    }
    fprintf(pFile, "%s\n", "");
    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "            </Cells>");
}

void VtuWriter::write_point_scalar_property(FILE *pFile, const std::vector<double> &data, std::string name)
{
    fprintf(pFile, "%s\n", "         <PointData Scalars=\"count\">");
    fprintf(pFile, "%s%s%s\n", "            <DataArray type=\"Float64\" Name=\"",
            name.c_str(),
            "\" "
            "format=\"ascii\">");
    for (size_t i = 0; i < size_t(data.size()); i++)
    {
        fprintf(pFile, "%e ", data[i]);
    }
    fprintf(pFile, "%s\n", "");

    fprintf(pFile, "%s\n", "            </DataArray>");

    fprintf(pFile, "%s\n", "         </PointData>");
}

void VtuWriter::write_cell_scalar_property(FILE *pFile, const std::vector<double> &data, std::string name)
{
    fprintf(pFile, "%s\n", "         <CellData Scalars=\"count\">");
    fprintf(pFile, "%s%s%s\n", "            <DataArray type=\"Float64\" Name=\"",
            name.c_str(),
            "\" "
            "format=\"ascii\">");
    for (size_t i = 0; i < size_t(data.size()); i++)
    {
        fprintf(pFile, "%e ", data[i]);
    }
    fprintf(pFile, "%s\n", "");

    fprintf(pFile, "%s\n", "            </DataArray>");

    fprintf(pFile, "%s\n", "         </CellData>");
}

void VtuWriter::write_footer(FILE *pFile)
{
    fprintf(pFile, "%s\n", "      </Piece>  ");
    fprintf(pFile, "%s\n", "   </UnstructuredGrid>");
    fprintf(pFile, "%s\n", "</VTKFile>");
}

bool VtuWriter::write_cells_to_vtu(std::string filename, const std::vector<double> &data)
{
    FILE *pFile = fopen(filename.c_str(), "w");
    write_header(pFile);
    write_vertices(pFile);
    write_cells(pFile);
    write_cell_scalar_property(pFile, data, "solution");
    write_footer(pFile);
    fclose(pFile);
    return true;
}

void VtuWriter::write_vector_field_property(FILE *pFile, const std::vector<Eigen::Vector2d> &vertex_values, std::string name)
{
    fprintf(pFile, "%s\n", "         <PointData Vectors=\"vector_field\">");
    fprintf(pFile, "%s%s%s\n", "            <DataArray type=\"Float64\" Name=\"",
            name.c_str(),
            "\" "
            "NumberOfComponents=\"3\" "
            "format=\"ascii\">");

    for (size_t i = 0; i < vertex_values.size(); i++)
    {
        Eigen::Vector2d vec = vertex_values[i];
        fprintf(pFile, "%e %e %e ", vec.x(), vec.y(), 0.0); // Assuming z-component is 0
    }

    fprintf(pFile, "%s\n", "");
    fprintf(pFile, "%s\n", "            </DataArray>");
    fprintf(pFile, "%s\n", "         </PointData>");
}

bool VtuWriter::write_to_vtu(std::string filename, const std::vector<double> &data, bool dimen)
{
    FILE *pFile = fopen(filename.c_str(), "w");
    write_header(pFile);
    if (dimen)
    {
        write_vertices(pFile, data);
    }
    else
    {
        write_vertices(pFile);
    }
    write_point_scalar_property(pFile, data, "enriched");
    write_cells(pFile);
    write_footer(pFile);
    fclose(pFile);
    return true;
}

bool VtuWriter::write_to_vtu(std::string filename, const std::vector<Eigen::Vector2d> &vertex_values, bool dimen)
{
    FILE *pFile = fopen(filename.c_str(), "w");
    write_header(pFile);
    write_vertices(pFile);
    // write_point_vector_property(pFile, vertex_values, "solution");
    write_vector_field_property(pFile, vertex_values, "enriched");
    write_cells(pFile);
    write_footer(pFile);
    fclose(pFile);
    return true;
}

bool VtuWriter::write_to_vtu(std::string filename)
{
    FILE *pFile = fopen(filename.c_str(), "w");
    write_header(pFile);
    write_vertices(pFile);
    write_cells(pFile);

    write_footer(pFile);
    fclose(pFile);
    return true;
}
