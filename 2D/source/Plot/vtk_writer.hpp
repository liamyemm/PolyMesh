#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <string>
#include <Eigen/Dense>
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkDoubleArray.h>
#include <vtkXMLPolyDataWriter.h>

#include "Mesh.hpp" // Include your Mesh class definition

class vtk_writer
{
public:
    vtk_writer(PolyMesh2D::CurvedMesh::Mesh *mesh, const std::vector<Eigen::Vector2d>& vertex_values);
    void write_to_vtk(const std::string& file_name);

private:
    PolyMesh2D::CurvedMesh::Mesh *mesh_;
    std::vector<Eigen::Vector2d> vertex_values_;
    vtkSmartPointer<vtkPoints> vtk_points_;
    vtkSmartPointer<vtkCellArray> cells_;
    vtkSmartPointer<vtkDoubleArray> vector_values_;
    vtkSmartPointer<vtkPolyData> poly_data_;

    void convert_mesh_to_vtk();
    void set_vertex_values();
};

#endif // VTK_WRITER_H
