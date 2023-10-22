#include "mesh_cutter.hpp"

#include <iostream>

using namespace PolyMesh2D;

MeshCutter::MeshCutter(CurvedMesh::Mesh *mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param, bool keep_internal, bool keep_external) : m_mesh(mesh), m_level_set(level_set), m_param(param), m_keep_internal(keep_internal), m_keep_external(keep_external)
{
    assert(m_keep_internal || m_keep_external);
}

void MeshCutter::cut_mesh()
{
    make_uncut_cells(m_mesh);
    make_cut_cells(m_mesh);
    make_boundary(m_mesh);
}

void MeshCutter::make_uncut_cells(CurvedMesh::Mesh *c_mesh)
{
    for (auto &cell : c_mesh->get_cells())
    {
        bool cell_internal = true;

        for (auto &vertex : cell->get_vertices())
        {
            if (m_level_set.value(vertex->coords()) < 1E-2 * m_mesh->h_max())
            {
                cell_internal = false;
                break;
            }
        }

        if (cell_internal)
        {
            if (!m_keep_internal)
            {
                c_mesh->remove_cell(cell);
            }
            continue;
        }

        bool cell_external = true;

        for (auto &vertex : cell->get_vertices())
        {
            if (m_level_set.value(vertex->coords()) > -1E-2 * m_mesh->h_max())
            {
                cell_external = false;
                break;
            }
        }

        // for external cells, we also check edge centers for extra robustness (although still not fully robust)
        for (auto &edge : cell->get_edges())
        {
            if (m_level_set.value(edge->center_mass()) > -1E-2 * m_mesh->h_max())
            // if (m_level_set.value(vertex->coords()) >= 0.0)
            {
                cell_external = false;
                break;
            }
        }

        if (cell_external)
        {
            if (!m_keep_external)
            {
                c_mesh->remove_cell(cell);
            }
            continue;
        }

        if (!cell_external && !cell_internal)
        {
            c_mesh->remove_cell(cell);
        }
    }

    for (auto &e : c_mesh->get_edges())
    {
        if (e->n_cells() == 0)
        {
            c_mesh->remove_edge(e);
        }
    }

    for (auto &v : c_mesh->get_vertices())
    {
        if(v->n_edges() == 0)
        {
            c_mesh->remove_vertex(v);
        }
    }
}

void MeshCutter::make_cut_cells(CurvedMesh::Mesh *c_mesh)
{
    std::vector<CurvedMesh::Edge *> current_edges = c_mesh->get_edges();

    // need to store curved edges and cells for later if curved edges overlap    
    std::vector<CurvedMesh::Edge *> curved_edges;
    std::vector<CurvedMesh::Cell *> curved_cells;

    for (auto &bdry_edge : current_edges)
    {
        if (bdry_edge->n_cells() != 1)
        {
            continue;
        }
        if (bdry_edge->is_boundary())
        {
            continue;
        }

        std::vector<CurvedMesh::Vertex *> new_cell_vertices;
        std::vector<CurvedMesh::Edge *> new_cell_edges;

        // allocate 4 edges and 4 vertices for the new cell
        new_cell_vertices.reserve(4); 
        new_cell_edges.reserve(4);

        new_cell_edges.push_back(bdry_edge);

        std::vector<double> t_vals; // parameter values of curved edge. Stored for later
        t_vals.reserve(2);

        for(auto & v : bdry_edge->get_vertices())
        {
            CurvedMesh::Vertex *v_new = nullptr;
            CurvedMesh::Edge *e_new = nullptr;

            double t = Functional::min_distance(m_param, v->coords());

            Eigen::Vector2d v_new_coord = m_param.value(t);
            t_vals.push_back(t);

            for (auto &vert : c_mesh->get_vertices()) 
            {
                if ((v_new_coord - vert->coords()).norm() < 1E-15)
                {
                    v_new = vert;
                    break;
                }
            }

            if(v_new != nullptr)
            {
                // v_new exist - check if e_new exists

                for (auto &edge : v->get_edges())
                {
                    if ((edge->vertex(0)->global_index() == v_new->global_index()) || (edge->vertex(1)->global_index() == v_new->global_index()))
                    {
                        e_new = edge;
                        break;
                    }
                }
            }
            else
            {
                // v_new does not exist - so make it

                v_new = new CurvedMesh::Vertex(c_mesh->n_vertices(), v_new_coord);
                c_mesh->add_vertex(v_new);
            }

            if(e_new == nullptr)
            {
                // e_new does not exist - so make it

                Functional::Curve edge_param = Functional::StraightLine(v->coords(), v_new_coord);
                e_new = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
                e_new->set_straight();
                c_mesh->add_edge(e_new);

                // add connections between new edges and vertices

                e_new->add_vertex(v);
                e_new->add_vertex(v_new);

                v->add_edge(e_new);
                v_new->add_edge(e_new);

                v->add_vertex(v_new);
                v_new->add_vertex(v);
            }

            // add to vectors of new edges and verts

            new_cell_edges.push_back(e_new);
            new_cell_vertices.push_back(v_new);
        }

        assert(new_cell_vertices.size() == 2);
        assert(new_cell_edges.size() == 3);
        assert(t_vals.size() == 2);

        Functional::Curve edge_param = Functional::restriction(m_param, t_vals[0], t_vals[1]);
        CurvedMesh::Edge *curved_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

        new_cell_edges.push_back(curved_edge);
        c_mesh->add_edge(curved_edge);

        for(auto & v : new_cell_vertices)
        {
            curved_edge->add_vertex(v);
            v->add_edge(curved_edge);
        }

        new_cell_vertices[0]->add_vertex(new_cell_vertices[1]);
        new_cell_vertices[1]->add_vertex(new_cell_vertices[0]);

        new_cell_vertices.push_back(bdry_edge->vertex(1));
        new_cell_vertices.push_back(bdry_edge->vertex(0));

        assert(new_cell_vertices.size() == 4);
        assert(new_cell_edges.size() == 4);

        CurvedMesh::Cell *curved_cell = new CurvedMesh::Cell(c_mesh->n_cells(), new_cell_edges);

        c_mesh->add_cell(curved_cell);

        for (auto &edge : new_cell_edges)
        {
            edge->add_cell(curved_cell);
        }
        for (auto &vert : new_cell_vertices)
        {
            vert->add_cell(curved_cell);
            curved_cell->add_vertex(vert);
        }

        curved_edges.push_back(curved_edge);
        curved_cells.push_back(curved_cell);
    }

    if(m_keep_internal && m_keep_external)
    {
        // curved edges overlap. Delete them and create non-overlapping ones.

        std::vector<CurvedMesh::Vertex *> verts;
        std::vector<double> t_vals;

        for (size_t iE = 0; iE < curved_edges.size(); ++iE)
        {
            if (std::find(verts.begin(), verts.end(), curved_edges[iE]->vertex(0)) == verts.end())
            {
                verts.push_back(curved_edges[iE]->vertex(0));
            }
            if (std::find(verts.begin(), verts.end(), curved_edges[iE]->vertex(1)) == verts.end())
            {
                verts.push_back(curved_edges[iE]->vertex(1));
            }
            c_mesh->remove_edge(curved_edges[iE]);
        }
        curved_edges.clear();

        for(auto & v : verts)
        {
            t_vals.push_back(Functional::min_distance(m_param, v->coords()));
        }

        assert(t_vals.size() == verts.size());

        for(size_t i = 0; i < t_vals.size(); ++i)
        {
            assert((m_param.value(t_vals[i]) - verts[i]->coords()).norm() < 1E-14);
        }

        for (size_t i = 0; i < t_vals.size(); ++i)
        {
            for (size_t j = i + 1; j < t_vals.size(); ++j)
            {
                if (t_vals[i] > t_vals[j])
                {
                    std::swap(t_vals[i], t_vals[j]);
                    std::swap(verts[i], verts[j]);
                }
            }
        }

        for (size_t i = 0; i < t_vals.size(); ++i)
        {
            size_t i_plus = (i + 1) % t_vals.size();

            if(i_plus == 0)
            {
                t_vals[i_plus] += m_param.tmax;
            }

            Functional::Curve edge_param = Functional::restriction(m_param, t_vals[i], t_vals[i_plus]);
            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            c_mesh->add_edge(new_edge);

            new_edge->add_vertex(verts[i]);
            new_edge->add_vertex(verts[i_plus]);

            verts[i]->add_edge(new_edge);
            verts[i_plus]->add_edge(new_edge);

            verts[i]->add_vertex(verts[i_plus]);
            verts[i_plus]->add_vertex(verts[i]);

            std::vector<CurvedMesh::Cell *> joint_cells;
            for (auto &cell : curved_cells)
            {
                std::vector<double> ts;
                for (auto &vert : cell->get_vertices())
                {
                    if (std::abs(m_level_set.value(vert->coords())) < 1E-15)
                    {
                        double t = Functional::min_distance(m_param, vert->coords());
                        ts.push_back(t);
                        assert((m_param.value(t) - vert->coords()).norm() < 1E-14);
                    }
                }
                double tmin = ts[0];
                double tmax = ts[0];

                size_t min_index = 0;
                for (size_t it = 1; it < ts.size(); ++it)
                {
                    double t = ts[it];
                    if (t < tmin)
                    {
                        tmin = t;
                        min_index = it;
                    }
                    if (t > tmax)
                    {
                        tmax = t;
                    }
                }

                bool tmax_fourth_quad = (tmax > m_param.tmin + 0.75 * (m_param.tmax - m_param.tmin));
                bool tmin_first_quad = (tmin < m_param.tmin + 0.25 * (m_param.tmax - m_param.tmin));
                if (tmin_first_quad && tmax_fourth_quad)
                {
                    ts[min_index] += m_param.tmax; // assumes periodic. will cause issues otherwise
                    tmax = ts[min_index]; // reset max to shifted min
                    tmin = ts[0];  // need to refind min
                    for (size_t it = 1; it < ts.size(); ++it)
                    {
                        double t = ts[it];
                        if (t < tmin)
                        {
                            tmin = t;
                        }
                    }
                }

                if (((t_vals[i] > tmin - 1E-15) && (t_vals[i] < tmax + 1E-15)) && ((t_vals[i_plus] > tmin - 1E-15) && (t_vals[i_plus] < tmax + 1E-15)))
                {
                    joint_cells.push_back(cell);
                }
            }

            assert(joint_cells.size() == 2); // should find exactly two cells

            for (auto &cell : joint_cells)
            {
                cell->add_edge(new_edge);
                new_edge->add_cell(cell);
                cell->set_edge_orientations();

                std::vector<CurvedMesh::Vertex *> cell_verts = cell->get_vertices();

                if (std::find(cell_verts.begin(), cell_verts.end(), verts[i]) == cell_verts.end())
                {
                    cell->add_vertex(verts[i]);
                    verts[i]->add_cell(cell);
                }

                if (std::find(cell_verts.begin(), cell_verts.end(), verts[i_plus]) == cell_verts.end())
                {
                    cell->add_vertex(verts[i_plus]);
                    verts[i_plus]->add_cell(cell);
                }
            }
        }
    }
}

void MeshCutter::make_boundary(CurvedMesh::Mesh *c_mesh)
{
    c_mesh->reset_boundary(); // empty all boundary and internal vectors, so we can reset them without duplicating.

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
            edge->set_boundary(false);
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
}
