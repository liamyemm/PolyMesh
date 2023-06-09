#include "mesh_convert.hpp"

namespace PolyMesh2D
{
    namespace MeshTransform
    {
        std::unique_ptr<CurvedMesh::Mesh> mesh_convert(const StraightMesh::Mesh *s_mesh)
        {
            std::unique_ptr<CurvedMesh::Mesh> c_mesh = std::make_unique<CurvedMesh::Mesh>();
            for (auto &cell : s_mesh->get_cells())
            {
                std::vector<CurvedMesh::Vertex *> new_vertices;

                // create vertices
                for (auto &cell_vert : cell->get_vertices())
                {
                    bool vert_exists = false;
                    for (auto &curved_vert : c_mesh->get_vertices())
                    {
                        if (cell_vert->coords() == curved_vert->coords()) // numerical precision comparison
                        {
                            vert_exists = true;
                            new_vertices.push_back(curved_vert);
                            break;
                        }
                    }
                    if (!vert_exists)
                    {
                        CurvedMesh::Vertex *new_vert = new CurvedMesh::Vertex(c_mesh->n_vertices(), cell_vert->coords());
                        new_vertices.push_back(new_vert);
                        c_mesh->add_vertex(new_vert);
                    }
                }

                std::vector<CurvedMesh::Edge *> new_edges;

                for (std::size_t i = 0; i < new_vertices.size(); ++i)
                {
                    std::size_t i_next = (i + 1) % new_vertices.size();

                    bool edge_exists = false;
                    if (new_vertices[i]->n_vertices() > 0)
                    {
                        std::vector<CurvedMesh::Vertex *> vlist = new_vertices[i]->get_vertices();
                        for (size_t j = 0; j < vlist.size(); j++)
                        {
                            if (vlist[j]->global_index() == new_vertices[i_next]->global_index()) // The edge exists
                            {
                                new_edges.push_back(new_vertices[i]->edge(j));
                                edge_exists = true;
                                break;
                            }
                        }
                    }
                    if (!edge_exists)
                    {
                        // edge does not exists - so create it
                        Functional::Curve edge_param = Functional::StraightLine(new_vertices[i]->coords(), new_vertices[i_next]->coords());
                        CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

                        new_edge->set_straight();

                        c_mesh->add_edge(new_edge);
                        new_edges.push_back(new_edge);

                        // add vertices to edge
                        new_edge->add_vertex(new_vertices[i]);
                        new_edge->add_vertex(new_vertices[i_next]);

                        // add edge to vertices
                        new_vertices[i]->add_edge(new_edge);
                        new_vertices[i_next]->add_edge(new_edge);

                        // add vertices to vertices
                        new_vertices[i]->add_vertex(new_vertices[i_next]);
                        new_vertices[i_next]->add_vertex(new_vertices[i]);
                    }
                }

                // create cell
                CurvedMesh::Cell *new_cell = new CurvedMesh::Cell(c_mesh->n_cells(), new_edges);
                c_mesh->add_cell(new_cell);

                for (auto &edge : new_edges)
                {
                    edge->add_cell(new_cell);
                }
                for (auto &vert : new_vertices)
                {
                    vert->add_cell(new_cell);
                    new_cell->add_vertex(vert);
                }
            }

            for (auto &edge : c_mesh->get_edges())
            {
                if (edge->n_cells() == 1) // If edge has only one cell, it is a boundary edge
                {
                    edge->set_boundary(true);          // set edge boundary to true
                    edge->cell(0)->set_boundary(true); // set cell boundary to true
                    edge->vertex(0)->set_boundary(true);
                    edge->vertex(1)->set_boundary(true);

                    c_mesh->add_b_edge(edge); // add edge to list of boundary edges
                }
                else
                {
                    c_mesh->add_i_edge(edge); // add edge to list of internal edges
                }
            }

            for (auto &cell : c_mesh->get_cells())
            {
                if (cell->is_boundary())
                {
                    c_mesh->add_b_cell(cell);
                }
                else
                {
                    c_mesh->add_i_cell(cell);
                }
            }

            for (auto &vert : c_mesh->get_vertices())
            {
                if (vert->is_boundary())
                {
                    c_mesh->add_b_vertex(vert);
                }
                else
                {
                    c_mesh->add_i_vertex(vert);
                }
            }

            return c_mesh;
        }

    }
}
