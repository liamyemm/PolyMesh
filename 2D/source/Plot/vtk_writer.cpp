#include "vtk_writer.hpp"
#include <array>
#include <functional>

vtk_writer::vtk_writer(PolyMesh2D::CurvedMesh::Mesh *mesh, const std::vector<Eigen::Vector2d>& vertex_values)
    : mesh_(mesh), vertex_values_(vertex_values)
{
    vtk_points_ = vtkSmartPointer<vtkPoints>::New();
    cells_ = vtkSmartPointer<vtkCellArray>::New();
    vector_values_ = vtkSmartPointer<vtkDoubleArray>::New();
    poly_data_ = vtkSmartPointer<vtkPolyData>::New();

    convert_mesh_to_vtk();
    // set_vertex_values();
}

void vtk_writer::convert_mesh_to_vtk()
{
    for (auto & vertex : mesh_->get_vertices())
    {
        double vtk_point[3] = {vertex->coords().x(), vertex->coords().y(), 0.0};
        vtk_points_->InsertNextPoint(vtk_point);
    }

    for (auto & cell : mesh_->get_cells())
    {
        vtkIdType cell_ids[cell->n_vertices()];
        for (size_t i = 0; i < cell->n_vertices(); ++i)
        {
            cell_ids[i] = cell->vertex(i)->global_index();
        }
        cells_->InsertNextCell(cell->n_vertices(), cell_ids);
    }
}

void vtk_writer::set_vertex_values()
{
    for (const Eigen::Vector2d& value : vertex_values_)
    {
        vector_values_->InsertNextTuple2(value.x(), value.y());
    }
}

void vtk_writer::write_to_vtk(const std::string& file_name)
{
    poly_data_->SetPoints(vtk_points_);
    poly_data_->SetPolys(cells_);
    // poly_data_->GetPointData()->SetVectors(vector_values_);

    vtkSmartPointer<vtkPolyData> dataset = vtkSmartPointer<vtkPolyData>::New();    
    dataset->ShallowCopy(poly_data_);

    vtkSmartPointer<vtkXMLDataSetWriter> writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(dataset);
    writer->Write();
}