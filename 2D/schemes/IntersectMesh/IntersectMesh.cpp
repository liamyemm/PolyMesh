#include "IntersectMesh.hpp"
#include <set>
#include <iterator>

int orientation(Eigen::Vector2d p, Eigen::Vector2d q, Eigen::Vector2d r)
{
    double val = (q(0) - p(0)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(0) - q(0));

    return Math::sgn(val);
}

using namespace PolyMesh2D;

MeshCutter::MeshCutter(StraightMesh::Mesh *s_mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param) : m_s_mesh(s_mesh), m_level_set(level_set), m_param(param)
{
}

std::unique_ptr<CurvedMesh::Mesh> MeshCutter::cut_mesh()
{
    std::unique_ptr<CurvedMesh::Mesh> c_mesh = std::make_unique<CurvedMesh::Mesh>();
    make_internal_cells(c_mesh.get());
    bool flag = true;
    if (c_mesh->n_cells() == 2)
    {
        flag = false;
        make_convex(c_mesh.get());
    }
    // make_convex(c_mesh.get());
    make_curved_cells(c_mesh.get());
    // if (flag)
    // {
    //     make_isotropic(c_mesh.get());
    // }
    // make_homogeneous(c_mesh.get());
    make_boundary(c_mesh.get());

    return c_mesh;
}

void MeshCutter::make_internal_cells(CurvedMesh::Mesh *c_mesh)
{
    for (auto &cell : m_s_mesh->get_cells())
    {
        bool cell_internal = true;

        for (auto &vertex : cell->get_vertices())
        {
            if (m_level_set.value(vertex->coords()) < 1E-1 * m_s_mesh->h_max())
            // if (m_level_set.value(vertex->coords()) <= 0.0)
            {
                cell_internal = false;
                break;
            }
        }

        if (cell_internal)
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
    }
}

void MeshCutter::make_convex(CurvedMesh::Mesh *c_mesh)
{
    std::vector<CurvedMesh::Edge *> non_convex_edges = c_mesh->get_edges();

    for (size_t i = 0; i < non_convex_edges.size(); ++i)
    {
        CurvedMesh::Edge *bdry_edge = non_convex_edges[i];
        // for (auto &bdry_edge : c_mesh->get_edges())
        // {
        if (bdry_edge->n_cells() != 1)
        {
            non_convex_edges.erase(non_convex_edges.begin() + i);
            --i;
            continue;
        }

        CurvedMesh::Vertex *v0 = bdry_edge->vertex(0);
        CurvedMesh::Vertex *v1 = bdry_edge->vertex(1);

        CurvedMesh::Vertex *v0_attached;
        CurvedMesh::Vertex *v1_attached;

        CurvedMesh::Edge *e0_attached;
        CurvedMesh::Edge *e1_attached;
        for (auto &b_edge : v0->get_edges())
        {
            if ((b_edge->global_index() != bdry_edge->global_index()) && (b_edge->n_cells() == 1))
            {
                v0_attached = b_edge->vertex(0)->global_index() == v0->global_index() ? b_edge->vertex(1) : b_edge->vertex(0);
                e0_attached = b_edge;
            }
        }

        for (auto &b_edge : v1->get_edges())
        {
            if ((b_edge->global_index() != bdry_edge->global_index()) && (b_edge->n_cells() == 1))
            {
                v1_attached = b_edge->vertex(0)->global_index() == v1->global_index() ? b_edge->vertex(1) : b_edge->vertex(0);
                e1_attached = b_edge;
            }
        }

        assert(v0_attached != nullptr);
        assert(v1_attached != nullptr);

        assert(e0_attached != nullptr);
        assert(e1_attached != nullptr);

        int o1 = orientation(v0->coords(), v1->coords(), v0_attached->coords());
        int o2 = orientation(v0->coords(), v1->coords(), v1_attached->coords());

        if (o1 != o2)
        {
            int good_orientation = orientation(v0->coords(), v1->coords(), bdry_edge->cell(0)->center_mass());

            assert(!((o1 == good_orientation) && (o2 == good_orientation)));

            double angle_1 = std::acos(((v0->coords() - v1->coords()).normalized()).dot((v0_attached->coords() - v0->coords()).normalized()));
            if (((v0->coords() - v1->coords()).normalized()).dot((v0_attached->coords() - v0->coords()).normalized()) >= 0.0)
            {
                angle_1 = Math::PI - angle_1;
            }
            double angle_2 = std::acos(((v1->coords() - v0->coords()).normalized()).dot((v1_attached->coords() - v1->coords()).normalized()));
            if (((v1->coords() - v0->coords()).normalized()).dot((v1_attached->coords() - v1->coords()).normalized()) >= 0.0)
            {
                angle_2 = Math::PI - angle_2;
            }

            bool ignore_v0_attached = (angle_1 > 0.74 * Math::PI);
            bool ignore_v1_attached = (angle_2 > 0.74 * Math::PI);

            if (ignore_v0_attached && ignore_v1_attached)
            {
                continue;
            }

            bool create_v1_v0_attached = !(o1 == good_orientation || ignore_v0_attached);
            bool create_v0_v1_attached = !(o2 == good_orientation || ignore_v1_attached);

            // assert(!(create_v1_v0_attached && create_v0_v1_attached));

            std::array<CurvedMesh::Vertex *, 2> edge_verts;
            std::vector<CurvedMesh::Vertex *> cell_verts({v0, v1});
            CurvedMesh::Edge *edge_attached;

            if (create_v1_v0_attached)
            {
                edge_verts[0] = v1;
                edge_verts[1] = v0_attached;

                edge_attached = e0_attached;
                cell_verts.push_back(v0_attached);
            }
            else if (create_v0_v1_attached)
            {
                edge_verts[0] = v0;
                edge_verts[1] = v1_attached;

                edge_attached = e1_attached;
                cell_verts.push_back(v1_attached);
            }
            else
            {
                non_convex_edges.erase(non_convex_edges.begin() + i);
                --i;
                continue;
            }

            Functional::Curve edge_param = Functional::StraightLine(edge_verts[0]->coords(), edge_verts[1]->coords());
            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            new_edge->set_straight();

            for (auto &vert : edge_verts)
            {
                vert->add_edge(new_edge);
                new_edge->add_vertex(vert);
            }

            edge_verts[0]->add_vertex(edge_verts[1]);
            edge_verts[1]->add_vertex(edge_verts[0]);

            c_mesh->add_edge(new_edge);

            std::vector<CurvedMesh::Edge *> cell_edges({new_edge, bdry_edge, edge_attached});

            CurvedMesh::Cell *new_cell = new CurvedMesh::Cell(c_mesh->n_cells(), cell_edges);

            c_mesh->add_cell(new_cell);

            for (auto &edge : cell_edges)
            {
                edge->add_cell(new_cell);
            }
            for (auto &vert : cell_verts)
            {
                vert->add_cell(new_cell);
                new_cell->add_vertex(vert);
            }
            i = 0;
            // goto convexify;
        }
        else
        {
            non_convex_edges.erase(non_convex_edges.begin() + i);
            --i;
        }
    }
}

void MeshCutter::make_curved_cells(CurvedMesh::Mesh *c_mesh)
{
    std::vector<CurvedMesh::Edge *> current_edges = c_mesh->get_edges();
    for (auto &bdry_edge : current_edges)
    {
        if (bdry_edge->n_cells() != 1)
        {
            continue;
        }

        // create two new vertices
        std::vector<CurvedMesh::Vertex *> new_vertices;
        std::vector<CurvedMesh::Edge *> new_edges;

        double t0 = Functional::min_distance(m_param, bdry_edge->vertex(0)->coords());
        Eigen::Vector2d new_vert_1_coord = m_param.value(t0);

        bool vert_1_exists = false;

        for (auto &b_vert : bdry_edge->vertex(0)->get_vertices())
        {
            if ((new_vert_1_coord - b_vert->coords()).norm() < 1E-14)
            {
                new_vertices.push_back(b_vert);
                vert_1_exists = true;
                break;
            }
        }

        if (!vert_1_exists)
        {
            // create vert and edge

            CurvedMesh::Vertex *new_vert = new CurvedMesh::Vertex(c_mesh->n_vertices(), new_vert_1_coord);
            new_vertices.push_back(new_vert);
            c_mesh->add_vertex(new_vert);

            Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(0)->coords(), new_vert_1_coord);
            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            new_edge->set_straight();

            new_edges.push_back(new_edge);
            c_mesh->add_edge(new_edge);

            new_edge->add_vertex(bdry_edge->vertex(0));
            new_edge->add_vertex(new_vert);

            bdry_edge->vertex(0)->add_edge(new_edge);
            new_vert->add_edge(new_edge);

            bdry_edge->vertex(0)->add_vertex(new_vert);
            new_vert->add_vertex(bdry_edge->vertex(0));
        }
        else
        {
            // edge exists
            for (auto &edge : bdry_edge->vertex(0)->get_edges())
            {
                if ((edge->vertex(0)->global_index() == new_vertices[0]->global_index()) || (edge->vertex(1)->global_index() == new_vertices[0]->global_index()))
                {
                    new_edges.push_back(edge);
                    break;
                }
            }
            assert(new_edges.size() == 1);
        }

        double t1 = Functional::min_distance(m_param, bdry_edge->vertex(1)->coords());
        Eigen::Vector2d new_vert_2_coord = m_param.value(t1);

        bool vert_2_exists = false;

        for (auto &b_vert : bdry_edge->vertex(1)->get_vertices())
        {
            if ((new_vert_2_coord - b_vert->coords()).norm() < 1E-14)
            {
                new_vertices.push_back(b_vert);
                vert_2_exists = true;
                break;
            }
        }

        if (!vert_2_exists)
        {
            // create vert and edge

            CurvedMesh::Vertex *new_vert = new CurvedMesh::Vertex(c_mesh->n_vertices(), new_vert_2_coord);
            new_vertices.push_back(new_vert);
            c_mesh->add_vertex(new_vert);

            Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(1)->coords(), new_vert_2_coord);
            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            new_edge->set_straight();

            new_edges.push_back(new_edge);
            c_mesh->add_edge(new_edge);

            new_edge->add_vertex(bdry_edge->vertex(1));
            new_edge->add_vertex(new_vert);

            bdry_edge->vertex(1)->add_edge(new_edge);
            new_vert->add_edge(new_edge);

            bdry_edge->vertex(1)->add_vertex(new_vert);
            new_vert->add_vertex(bdry_edge->vertex(1));
        }
        else
        {
            // edge exists
            for (auto &edge : bdry_edge->vertex(1)->get_edges())
            {
                if ((edge->vertex(0)->global_index() == new_vertices[1]->global_index()) || (edge->vertex(1)->global_index() == new_vertices[1]->global_index()))
                {
                    new_edges.push_back(edge);
                    break;
                }
            }
            assert(new_edges.size() == 2);
        }

        Functional::Curve edge_param = Functional::restriction(m_param, t0, t1);
        CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

        // Functional::Curve edge_param = Functional::StraightLine(new_vert_1_coord, new_vert_2_coord);
        // CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
        // new_edge->set_straight();

        new_edges.push_back(new_edge);
        c_mesh->add_edge(new_edge);

        new_edge->add_vertex(new_vertices[0]);
        new_edge->add_vertex(new_vertices[1]);

        new_vertices[0]->add_edge(new_edge);
        new_vertices[1]->add_edge(new_edge);

        new_vertices[0]->add_vertex(new_vertices[1]);
        new_vertices[1]->add_vertex(new_vertices[0]);

        new_vertices.push_back(bdry_edge->vertex(1));
        new_vertices.push_back(bdry_edge->vertex(0));

        new_edges.push_back(bdry_edge);
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
}

void MeshCutter::make_isotropic(CurvedMesh::Mesh *c_mesh)
{
    // std::vector<CurvedMesh::Cell *> mergable_cells = c_mesh->get_cells();

    // for (size_t i = 0; i < mergable_cells.size(); ++i)
    // {
    //     CurvedMesh::Cell *cell = mergable_cells[i];
    for (size_t i = 0; i < c_mesh->n_cells(); ++i)
    {
        CurvedMesh::Cell *cell = c_mesh->cell(i);
        double diam = cell->diam();
        double cell_measure = cell->measure();
        double bdry_measure = 0.0;
        for (auto &e : cell->get_edges())
        {
            bdry_measure += e->measure();
        }

        double eps = 0.075;
        double isotropy = cell_measure / bdry_measure;
        if (isotropy < eps * diam)
        {
            CurvedMesh::Edge *longest_edge = ( (cell->edge(0)->n_cells() == 1) ? cell->edge(1) : cell->edge(0));
            if (longest_edge->n_cells() == 1)
            {
                longest_edge = cell->edge(2);
            }
            // errors if cell has three or more boundary edges and these are the longest edges.
            for (size_t iE = 1; iE < cell->n_edges(); ++iE)
            {
                if (cell->edge(iE)->n_cells() == 1)
                {
                    continue;
                }
                if (cell->edge(iE)->measure() > longest_edge->measure())
                {
                    longest_edge = cell->edge(iE);
                }
            }
            assert(longest_edge->n_cells() == 2);
            CurvedMesh::Cell *merging_cell = (longest_edge->cell(0)->global_index() == cell->global_index() ? longest_edge->cell(1) : longest_edge->cell(0));

            int n_shared_edges = 0;
            for (auto &e : merging_cell->get_edges())
            {
                auto vec = cell->get_edges();
                if (std::find(vec.begin(), vec.end(), e) != vec.end())
                {
                    ++n_shared_edges;
                }
            }
            if (n_shared_edges != 1)
            {
                // std::cout << "\n" << n_shared_edges << "\n";
                continue;
            }

            int n_shared_verts = 0;
            for (auto &v : merging_cell->get_vertices())
            {
                auto vec = cell->get_vertices();
                if (std::find(vec.begin(), vec.end(), v) != vec.end())
                {
                    ++n_shared_verts;
                }
            }
            if (n_shared_verts != 2)
            {
                std::cout << "\n"
                          << n_shared_verts << "\n";
                continue;
            }

            std::vector<CurvedMesh::Vertex *> combined_verts = cell->get_vertices();
            std::vector<CurvedMesh::Edge *> combined_edges = cell->get_edges();

            std::vector<CurvedMesh::Vertex *> merging_cell_verts = merging_cell->get_vertices();
            std::vector<CurvedMesh::Edge *> merging_cell_edges = merging_cell->get_edges();

            combined_verts.insert(combined_verts.end(), merging_cell_verts.begin(), merging_cell_verts.end());
            combined_edges.insert(combined_edges.end(), merging_cell_edges.begin(), merging_cell_edges.end());

            std::set<CurvedMesh::Vertex *> set_verts(combined_verts.begin(), combined_verts.end()); // convert to set to erase duplicates
            combined_verts.assign(set_verts.begin(), set_verts.end());                              // convert back to vector

            std::set<CurvedMesh::Edge *> set_edges(combined_edges.begin(), combined_edges.end()); // convert to set to erase duplicates
            combined_edges.assign(set_edges.begin(), set_edges.end());                            // convert back to vector

            combined_edges.erase(std::find(combined_edges.begin(), combined_edges.end(), longest_edge)); // remove merging edge

            assert(combined_edges.size() == cell->n_edges() + merging_cell->n_edges() - 2);
            assert(combined_verts.size() == cell->n_vertices() + merging_cell->n_vertices() - 2);

            CurvedMesh::Cell *new_cell = new CurvedMesh::Cell(c_mesh->n_cells(), combined_edges);

            c_mesh->add_cell(new_cell);

            c_mesh->remove_cell(cell);
            c_mesh->remove_cell(merging_cell);
            c_mesh->remove_edge(longest_edge);

            // double new_cell_measure = new_cell->measure();
            // double new_bdry_measure = 0.0;
            // for (auto &e : new_cell->get_edges())
            // {
            //     new_bdry_measure += e->measure();
            // }

            // double new_isotropy = new_cell_measure / new_bdry_measure;

            // if(new_isotropy < isotropy)
            // {
            //     continue;
            // }

            // mergable_cells.erase(std::find(mergable_cells.begin(), mergable_cells.end(), cell));
            // mergable_cells.erase(std::find(mergable_cells.begin(), mergable_cells.end(), merging_cell));
            // i -= 2;
            // mergable_cells.push_back(new_cell);

            for (auto &edge : combined_edges)
            {
                edge->add_cell(new_cell);
            }
            for (auto &vert : combined_verts)
            {
                vert->add_cell(new_cell);
                new_cell->add_vertex(vert);
            }
            i = 0;
        }
        // else
        // {
        // mergable_cells.erase(mergable_cells.begin() + i);
        // --i;
        // }
    }
}

void MeshCutter::make_homogeneous(CurvedMesh::Mesh *c_mesh)
{
    std::vector<CurvedMesh::Cell *> cuttable_cells = c_mesh->get_cells();

    for (size_t i = 0; i < cuttable_cells.size(); ++i)
    {
        CurvedMesh::Cell *cell = cuttable_cells[i];

        double tol = 1.5;
        if (cell->diam() > tol * m_s_mesh->h_max())
        {
            CurvedMesh::Edge *longest_edge = cell->edge(0);
            for (size_t iE = 1; iE < cell->n_edges(); ++iE)
            {
                if (cell->edge(iE)->measure() > longest_edge->measure())
                {
                    longest_edge = cell->edge(iE);
                }
            }
            double tmin = longest_edge->parameterisation().tmin;
            double tmax = longest_edge->parameterisation().tmax;
            double t_half = 0.5 * (tmin + tmax);
            Eigen::Vector2d point_half_way = longest_edge->parameterisation().value(t_half);

            CurvedMesh::Vertex *vert_to_connect = cell->vertex(0);
            std::function<bool(size_t)> is_valid;

            if (cell->n_vertices() < 5)
            {
                // find vertex not attached to edge that is a maximum distance from point_half_way
                is_valid = [&longest_edge](const size_t id) -> bool
                {
                    bool valid = (id != longest_edge->vertex(0)->global_index());
                    valid = (valid && (id != longest_edge->vertex(1)->global_index()));
                    return valid;
                };
            }
            else
            {
                // find vertex not attached to edge and is not adjacent to either edge vert that is a maximum distance from point_half_way
                is_valid = [&longest_edge, &c_mesh](const size_t id) -> bool
                {
                    CurvedMesh::Vertex *v1 = longest_edge->vertex(0);
                    CurvedMesh::Vertex *v2 = longest_edge->vertex(1);
                    bool valid = (id != v1->global_index());
                    valid = (valid && (id != v2->global_index()));

                    if (valid)
                    {
                        std::vector<CurvedMesh::Vertex *> test_verts = c_mesh->vertex(id)->get_vertices();
                        if (std::find(test_verts.begin(), test_verts.end(), v1) != test_verts.end())
                        {
                            valid = false;
                        }
                        if (std::find(test_verts.begin(), test_verts.end(), v2) != test_verts.end())
                        {
                            valid = false;
                        }
                    }
                    return valid;
                };
            }

            for (auto &vert : cell->get_vertices())
            {
                if (!is_valid(vert->global_index()))
                {
                    continue;
                }
                if (!is_valid(vert_to_connect->global_index()))
                {
                    vert_to_connect = vert;
                    continue;
                }
                if ((point_half_way - vert->coords()).norm() > (point_half_way - vert_to_connect->coords()).norm())
                {
                    vert_to_connect = vert;
                }
            }

            assert(is_valid(vert_to_connect->global_index()));

            std::vector<CurvedMesh::Vertex *> ordered_verts = longest_edge->get_vertices();

            size_t pos = 0;
            while (ordered_verts.size() < cell->n_vertices())
            {
                if (std::find(ordered_verts.begin(), ordered_verts.end(), cell->vertex(pos)) == ordered_verts.end())
                {
                    // not in ordered verts
                    std::vector<CurvedMesh::Vertex *> test_verts = cell->vertex(pos)->get_vertices();
                    if (std::find(test_verts.begin(), test_verts.end(), ordered_verts[ordered_verts.size() - 1]) != test_verts.end())
                    {
                        // connects to final vert of ordered verts
                        ordered_verts.push_back(cell->vertex(pos));
                    }
                }
                pos = (pos + 1) % cell->n_vertices();
            }

            assert(ordered_verts.size() == cell->n_vertices());
            ordered_verts.erase(ordered_verts.begin() + 0);
            ordered_verts.push_back(longest_edge->vertex(0));

            auto iter = std::find(ordered_verts.begin(), ordered_verts.end(), vert_to_connect);

            assert(iter != ordered_verts.end());

            std::vector<CurvedMesh::Vertex *> vert_set_1(ordered_verts.begin(), iter + 1);
            std::vector<CurvedMesh::Vertex *> vert_set_2(iter, ordered_verts.end());

            assert(vert_set_1.size() + vert_set_2.size() == cell->n_vertices() + 1);

            CurvedMesh::Vertex *new_vert = new CurvedMesh::Vertex(c_mesh->n_vertices(), point_half_way);
            c_mesh->add_vertex(new_vert);

            new_vert->add_vertex(vert_to_connect);
            vert_to_connect->add_vertex(new_vert);

            new_vert->add_vertex(longest_edge->vertex(0));
            longest_edge->vertex(0)->add_vertex(new_vert);

            new_vert->add_vertex(longest_edge->vertex(1));
            longest_edge->vertex(1)->add_vertex(new_vert);

            vert_set_1.push_back(new_vert);
            vert_set_2.push_back(new_vert);

            Functional::Curve edge_cut_1_param;
            Functional::Curve edge_cut_2_param;

            if ((longest_edge->parameterisation().value(tmin) - longest_edge->vertex(0)->coords()).norm() < 1E-12)
            {
                assert((longest_edge->parameterisation().value(tmax) - longest_edge->vertex(1)->coords()).norm() < 1E-12);
                edge_cut_1_param = Functional::restriction(longest_edge->parameterisation(), tmin, t_half);
                edge_cut_2_param = Functional::restriction(longest_edge->parameterisation(), t_half, tmax);
            }
            else
            {
                assert((longest_edge->parameterisation().value(tmin) - longest_edge->vertex(1)->coords()).norm() < 1E-12);
                assert((longest_edge->parameterisation().value(tmax) - longest_edge->vertex(0)->coords()).norm() < 1E-12);
                edge_cut_2_param = Functional::restriction(longest_edge->parameterisation(), tmin, t_half);
                edge_cut_1_param = Functional::restriction(longest_edge->parameterisation(), t_half, tmax);
            }

            CurvedMesh::Edge *edge_cut_1 = new CurvedMesh::Edge(c_mesh->n_edges(), edge_cut_1_param);
            c_mesh->add_edge(edge_cut_1);

            CurvedMesh::Edge *edge_cut_2 = new CurvedMesh::Edge(c_mesh->n_edges(), edge_cut_2_param);
            c_mesh->add_edge(edge_cut_2);

            edge_cut_1->add_vertex(longest_edge->vertex(0));
            longest_edge->vertex(0)->add_edge(edge_cut_1);

            new_vert->add_edge(edge_cut_1);
            edge_cut_1->add_vertex(new_vert);

            edge_cut_2->add_vertex(longest_edge->vertex(1));
            longest_edge->vertex(1)->add_edge(edge_cut_2);

            new_vert->add_edge(edge_cut_2);
            edge_cut_2->add_vertex(new_vert);

            if (longest_edge->n_cells() == 2)
            {
                // internal edge
                // get non boundary cell
                CurvedMesh::Cell *other_cell = (cell->global_index() == longest_edge->cell(0)->global_index() ? longest_edge->cell(1) : longest_edge->cell(0));

                edge_cut_1->add_cell(other_cell);
                edge_cut_2->add_cell(other_cell);
                other_cell->add_edge(edge_cut_1);
                other_cell->add_edge(edge_cut_2);

                other_cell->add_vertex(new_vert);
                new_vert->add_cell(other_cell);
            }

            if (longest_edge->is_straight())
            {
                edge_cut_1->set_straight();
                edge_cut_2->set_straight();
            }

            std::vector<CurvedMesh::Edge *> all_edges = cell->get_edges();
            all_edges.push_back(edge_cut_1);
            all_edges.push_back(edge_cut_2);
            all_edges.erase(std::find(all_edges.begin(), all_edges.end(), longest_edge));

            c_mesh->remove_cell(cell);
            c_mesh->remove_edge(longest_edge);

            std::vector<CurvedMesh::Edge *> edge_set_1;
            std::vector<CurvedMesh::Edge *> edge_set_2;
            for (auto &e : all_edges)
            {
                bool found_e0_set_1 = (std::find(vert_set_1.begin(), vert_set_1.end(), e->vertex(0)) != vert_set_1.end());
                bool found_e0_set_2 = (std::find(vert_set_2.begin(), vert_set_2.end(), e->vertex(0)) != vert_set_2.end());
                bool found_e1_set_1 = (std::find(vert_set_1.begin(), vert_set_1.end(), e->vertex(1)) != vert_set_1.end());
                bool found_e1_set_2 = (std::find(vert_set_2.begin(), vert_set_2.end(), e->vertex(1)) != vert_set_2.end());

                if (found_e0_set_1 && found_e1_set_1)
                {
                    edge_set_1.push_back(e);
                }
                else if (found_e0_set_2 && found_e1_set_2)
                {
                    edge_set_2.push_back(e);
                }
            }

            Functional::Curve joint_edge_param = Functional::StraightLine(point_half_way, vert_to_connect->coords());
            CurvedMesh::Edge *joint_edge = new CurvedMesh::Edge(c_mesh->n_edges(), joint_edge_param);
            c_mesh->add_edge(joint_edge);

            joint_edge->set_straight();

            joint_edge->add_vertex(vert_to_connect);
            joint_edge->add_vertex(new_vert);
            vert_to_connect->add_edge(joint_edge);
            new_vert->add_edge(joint_edge);

            edge_set_1.push_back(joint_edge);
            edge_set_2.push_back(joint_edge);

            CurvedMesh::Cell *new_cell_1 = new CurvedMesh::Cell(c_mesh->n_cells(), edge_set_1);
            CurvedMesh::Cell *new_cell_2 = new CurvedMesh::Cell(c_mesh->n_cells(), edge_set_2);

            c_mesh->add_cell(new_cell_1);
            c_mesh->add_cell(new_cell_2);

            for (auto &edge : edge_set_1)
            {
                edge->add_cell(new_cell_1);
            }

            for (auto &edge : edge_set_2)
            {
                edge->add_cell(new_cell_2);
            }

            for (auto &vert : vert_set_1)
            {
                vert->add_cell(new_cell_1);
                new_cell_1->add_vertex(vert);
            }

            for (auto &vert : vert_set_2)
            {
                vert->add_cell(new_cell_2);
                new_cell_2->add_vertex(vert);
            }

            cuttable_cells.push_back(new_cell_1);
            cuttable_cells.push_back(new_cell_2);
            cuttable_cells.erase(cuttable_cells.begin() + i);
            --i;
        }
        else
        {
            cuttable_cells.erase(cuttable_cells.begin() + i);
            --i;
        }
    }
}

void MeshCutter::make_boundary(CurvedMesh::Mesh *c_mesh)
{
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
}
