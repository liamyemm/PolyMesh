#include "mesh_builder.hpp"

using namespace PolyMesh3D::StraightMesh;

MeshBuilder::MeshBuilder(const std::string mesh_file) : _mesh_file(mesh_file) {}

std::unique_ptr<Mesh> MeshBuilder::build_the_mesh()
{
    MeshReaderRF mesh_reader(_mesh_file);

    std::vector<std::vector<double>> vertices;
    std::vector<std::vector<std::vector<size_t>>> cells;

    mesh_reader.read_node_file(vertices);
    mesh_reader.read_ele_file(cells);

    if (vertices.size() > 0 && cells.size() > 0)
    {
        std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>(); // make a pointer to the mesh so that it outlives the builder

        std::cout << "     Mesh: ";

        // Create vertices
        for (auto &v : vertices)
        {
            VectorRd vert(v[0], v[1], v[2]);
            Vertex *vertex = new Vertex(mesh->n_vertices(), vert);
            mesh->add_vertex(vertex);
        }

        // Create cells
        double total_vol = 0.0;
        for (auto &c : cells)
        {
            // build edges and faces of cell
            std::vector<Face *> cell_faces;
            std::vector<Edge *> cell_edges;
            // std::vector<VectorRd> cell_vertex_coords;
            std::vector<std::size_t> cell_vertex_ids;
            for (auto &f : c)
            {
                std::vector<std::size_t> face_vertex_ids;
                // std::vector<VectorRd> face_vertex_coords;
                for (auto &vertID : f)
                {
                    face_vertex_ids.push_back(vertID);
                    VectorRd coord({vertices[vertID][0], vertices[vertID][1], vertices[vertID][2]});
                    // face_vertex_coords.push_back(coord);
                    if (std::find(cell_vertex_ids.begin(), cell_vertex_ids.end(), vertID) == cell_vertex_ids.end())
                    {
                        cell_vertex_ids.push_back(vertID);
                        // cell_vertex_coords.push_back(coord);
                    }
                }

                // make edges
                // bool make_face_flag = true;
                std::vector<Edge *> face_edges;
                // std::vector<std::size_t> need_to_add;
                // Face *face;
                bool face_does_not_exist = false;
                for (std::size_t i = 0; i < face_vertex_ids.size(); ++i)
                {
                    std::array<size_t, 2> e({face_vertex_ids[i], face_vertex_ids[(i + 1) % face_vertex_ids.size()]});

                    bool edge_not_found = true;

                    std::vector<Vertex *> vlist = mesh->vertex(e[0])->get_vertices();
                    for (size_t j = 0; j < vlist.size(); j++)
                    {
                        if (vlist[j]->global_index() == e[1]) // The edge exists in the mesh
                        {
                            Edge *edge = mesh->vertex(e[0])->edge(j);
                            // need_to_add.push_back(face_edges.size());
                            face_edges.push_back(edge);
                            if (std::find(cell_edges.begin(), cell_edges.end(), edge) == cell_edges.end())
                            {
                                cell_edges.push_back(edge);
                            }
                            edge_not_found = false;
                            break;
                        }
                    }

                    if (edge_not_found)
                    {
                        face_does_not_exist = true;
                        Vertex *vertex1 = mesh->vertex(e[0]);
                        Vertex *vertex2 = mesh->vertex(e[1]);

                        Edge *edge = new Edge(mesh->n_edges(), {vertex1, vertex2});

                        cell_edges.push_back(edge);
                        face_edges.push_back(edge);

                        mesh->add_edge(edge);

                        vertex1->add_edge(edge);
                        vertex2->add_edge(edge);

                        vertex1->add_vertex(vertex2);
                        vertex2->add_vertex(vertex1);
                    }
                }

                Face *face;
                bool face_exists = !face_does_not_exist;
                if(face_exists)
                {
                    // face might still not exist. Need to confirm.
                    face_exists = false;
                    for(auto & existing_face : mesh->get_faces())
                    {
                        auto existing_face_edges = existing_face->get_edges();
                        bool is_this_face = true;
                        for(auto & edge : face_edges)
                        {
                            if(std::find(existing_face_edges.begin(), existing_face_edges.end(), edge) == existing_face_edges.end())
                            {
                                is_this_face = false;
                                break;
                            }
                        }
                        if(is_this_face)
                        {
                            face_exists = true;
                            face = existing_face;
                            break;
                        }
                    }
                }

                if (!face_exists) 
                {
                    // make face
                    face = new Face(mesh->n_faces(), face_edges);
                    mesh->add_face(face);
                    for (auto &e : face_edges)
                    {
                        e->add_face(face);
                    }
                    for (auto &vertID : face_vertex_ids)
                    {
                        face->add_vertex(mesh->vertex(vertID));
                        mesh->vertex(vertID)->add_face(face);
                    }
                }
                cell_faces.push_back(face);
            }

            Cell *cell = new Cell(mesh->n_cells(), cell_faces);
            total_vol += cell->measure();
            mesh->add_cell(cell);

            for (auto &face : cell_faces)
            {
                face->add_cell(cell);
            }

            for (auto &edge : cell_edges)
            {
                cell->add_edge(edge);
                edge->add_cell(cell);
            }

            for (auto &vertID : cell_vertex_ids)
            {
                cell->add_vertex(mesh->vertex(vertID));
                mesh->vertex(vertID)->add_cell(cell);
            }
        }

        // // build boundary
        this->build_boundary(mesh.get());

        std::cout << "added " << mesh->n_cells() << " cells; Total volume = " << total_vol << std::endl;
        return mesh;
    }
    else
    {
        throw "     Cannot build mesh. Check input file\n";
    }
    return NULL;
}

void MeshBuilder::build_boundary(Mesh *mesh)
{
    // Here we fill in the _boundary variables of the cells and vertices, and the lists of boundary
    // edges, cells and vertices
    for (auto &face : mesh->get_faces())
    {
        std::vector<Cell *> cells = face->get_cells();
        if (cells.size() == 1)
        {
            // The cell has a boundary face, so it is a boundary cell
            cells[0]->set_boundary(true);
            face->set_boundary(true);
            // mesh->add_b_cell(cells[0]);
            mesh->add_b_face(face);

            std::vector<Vertex *> b_verts = mesh->get_b_vertices();
            for (auto &v : face->get_vertices())
            {
                if (std::find(b_verts.begin(), b_verts.end(), v) == b_verts.end()) // if vert not in b_verts, add to b_verts
                {
                    v->set_boundary(true);
                    mesh->add_b_vertex(v);
                }
            }

            std::vector<Edge *> b_edges = mesh->get_b_edges();
            for (auto &e : face->get_edges())
            {
                if (std::find(b_edges.begin(), b_edges.end(), e) == b_edges.end()) // if edge not in b_edges, add to b_edges
                {
                    e->set_boundary(true);
                    mesh->add_b_edge(e);
                }
            }
        }
        else
        {
            mesh->add_i_face(face);
        }
    }

    for (auto &cell : mesh->get_cells())
    {
        if (!(cell->is_boundary()))
        {
            mesh->add_i_cell(cell);
        }
        else
        {
            mesh->add_b_cell(cell);
        }
    }
    for (auto &edge : mesh->get_edges())
    {
        if (!(edge->is_boundary()))
        {
            mesh->add_i_edge(edge);
        }
    }
    for (auto &vertex : mesh->get_vertices())
    {
        if (!(vertex->is_boundary()))
        {
            mesh->add_i_vertex(vertex);
        }
    }
}
