#include "IntersectCurvedMesh.hpp"
#include <set>
#include <iterator>

int orient(Eigen::Vector2d p, Eigen::Vector2d q, Eigen::Vector2d r)
{
    double val = (q(0) - p(0)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(0) - q(0));

    return Math::sgn(val);
}

using namespace PolyMesh2D;

// MeshCutter2::MeshCutter2(StraightMesh::Mesh *s_mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param) : m_s_mesh(s_mesh), m_level_set(level_set), m_param(param)
// {
// }

MeshIntersect::MeshIntersect(CurvedMesh::Mesh *s_mesh, const Functional::ScalarFunction2D &level_set, const Functional::Curve &param) : m_s_mesh(s_mesh), m_level_set(level_set), m_param(param)
{
    no_external = true;
    make_straight = false;
}

// std::unique_ptr<CurvedMesh::Mesh> MeshCutter2::cut_mesh()
void MeshIntersect::cut_mesh()
{

    CurvedMesh::Mesh * c_mesh = m_s_mesh;
    // std::unique_ptr<CurvedMesh::Mesh> c_mesh = std::make_unique<CurvedMesh::Mesh>();
    make_internal_cells(c_mesh);
    // make_convex(c_mesh.get());
    make_curved_cells(c_mesh);
    // if(m_s_mesh->h_max() < 0.5)
    // {
    //     make_isotropic(c_mesh.get());
    // }
    // make_boundary(c_mesh);

    // return c_mesh;
}

void MeshIntersect::make_convex(CurvedMesh::Mesh *c_mesh)
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

        int o1 = orient(v0->coords(), v1->coords(), v0_attached->coords());
        int o2 = orient(v0->coords(), v1->coords(), v1_attached->coords());

        if (o1 != o2)
        {
            int good_orientation = orient(v0->coords(), v1->coords(), bdry_edge->cell(0)->center_mass());

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

void MeshIntersect::make_internal_cells(CurvedMesh::Mesh *c_mesh)
{
    for (auto &cell : c_mesh->get_cells())
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

        bool cell_external = true;

        for (auto &vertex : cell->get_vertices())
        {
            if (m_level_set.value(vertex->coords()) > -1E-1 * m_s_mesh->h_max())
            // if (m_level_set.value(vertex->coords()) >= 0.0)
            {
                cell_external = false;
                break;
            }
        }

        if (cell_external)
        {
            no_external = false;
        }

        if ( !(cell_internal || cell_external) )
        {
            c_mesh->remove_cell(cell);
        }
    }

    for (auto &e : c_mesh->get_edges())
    {
        if(e->n_cells() == 0)
        {
            c_mesh->remove_edge(e);
        }
    }

    // for (auto &v : c_mesh->get_vertices())
    // {
    //     if(v->n_edges() == 0)
    //     {
    //         c_mesh->remove_vertex(v);
    //     }
    // }
}

void MeshIntersect::make_curved_cells(CurvedMesh::Mesh *c_mesh)
{
    std::vector<CurvedMesh::Edge *> current_edges = c_mesh->get_edges();

    // std::vector<CurvedMesh::Vertex *> verts;
    // std::vector<double> t_vals;

    std::vector<CurvedMesh::Edge *> curved_edges;

    for (auto &bdry_edge : current_edges)
    {
        if (bdry_edge->n_cells() != 1)
        {
            continue;
        }

        // double cx = bdry_edge->center_mass()(0);
        // double cy = bdry_edge->center_mass()(1);

        // if (std::abs(cx - 0.0) < 1E-15 || std::abs(cx - 1.0) < 1E-15 || std::abs(cy - 0.0) < 1E-15 || std::abs(cy - 1.0) < 1E-15)
        // {
        //     continue;
        // }

        // if (std::abs(m_level_set.value(bdry_edge->vertex(0)->coords())) > 0.1 || std::abs(m_level_set.value(bdry_edge->vertex(1)->coords())) > 0.1)
        // {
        //     continue;
        // }

        if (bdry_edge->is_boundary())
        {
            continue;
        }

        // std::cout << "dafuq?\n";

        // if( m_level_set.value(bdry_edge->center_mass()) < 0.0 )
        // {
        //     continue;
        // }

        // create two new vertices
        std::vector<CurvedMesh::Vertex *> new_vertices;
        std::vector<CurvedMesh::Edge *> new_edges;

        double t0 = Functional::min_distance(m_param, bdry_edge->vertex(0)->coords());

        // double v0x = bdry_edge->vertex(0)->coords()(0);
        // double v0y = bdry_edge->vertex(0)->coords()(1);

        // if (std::abs(v0x) < 1E-15)
        // {
        //     t0 = Math::PI / 2.0;
        //     // t0 = Functional::min_distance(m_param, Eigen::Vector2d(0.0, std::sqrt(0.9375) - 0.25));
        // }
        // else if (std::abs(v0y) < 1E-15)
        // {
        //     t0 = 0.0;
        //     // t0 = Functional::min_distance(m_param, Eigen::Vector2d(std::sqrt(0.9375) - 0.25, 0.0));
        // }
        // else
        // {
        t0 = Functional::min_distance(m_param, bdry_edge->vertex(0)->coords());
        // }

        Eigen::Vector2d new_vert_1_coord = m_param.value(t0);

        bool vert_1_exists = false;

        // for (auto &b_vert : bdry_edge->vertex(0)->get_vertices()) // assumes no hanging nodes
        for (auto &b_vert : c_mesh->get_vertices())
        {
            if ((new_vert_1_coord - b_vert->coords()).norm() < 1E-15)
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

            c_mesh->add_i_vertex(new_vert);

            Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(0)->coords(), new_vert_1_coord);
            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            c_mesh->add_i_edge(new_edge);

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
            if (new_edges.size() != 1)
            {
                Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(0)->coords(), new_vert_1_coord);
                CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

                c_mesh->add_i_edge(new_edge);

                new_edge->set_straight();

                new_edges.push_back(new_edge);
                c_mesh->add_edge(new_edge);

                new_edge->add_vertex(bdry_edge->vertex(0));
                new_edge->add_vertex(new_vertices[0]);

                bdry_edge->vertex(0)->add_edge(new_edge);
                new_vertices[0]->add_edge(new_edge);

                bdry_edge->vertex(0)->add_vertex(new_vertices[0]);
                new_vertices[0]->add_vertex(bdry_edge->vertex(0));
            }
            assert(new_edges.size() == 1);
        }

        double t1 = Functional::min_distance(m_param, bdry_edge->vertex(1)->coords());

        // double v1x = bdry_edge->vertex(1)->coords()(0);
        // double v1y = bdry_edge->vertex(1)->coords()(1);

        // if (std::abs(v1x) < 1E-15)
        // {
        //     t1 = Math::PI / 2.0;
        //     // t1 = Functional::min_distance(m_param, Eigen::Vector2d(0.0, std::sqrt(0.9375) - 0.25));
        // }
        // else if (std::abs(v1y) < 1E-15)
        // {
        //     t1 = 0.0;
        //     // t1 = Functional::min_distance(m_param, Eigen::Vector2d(std::sqrt(0.9375) - 0.25, 0.0));
        // }
        // else
        // {
        // t1 = Functional::min_distance(m_param, bdry_edge->vertex(1)->coords());
        // }

        assert(std::abs(t0 - t1) > 1E-10);

        Eigen::Vector2d new_vert_2_coord = m_param.value(t1);

        // double t1 = Functional::min_distance(m_param, bdry_edge->vertex(1)->coords());
        // Eigen::Vector2d new_vert_2_coord = m_param.value(t1);

        bool vert_2_exists = false;

        // for (auto &b_vert : bdry_edge->vertex(1)->get_vertices()) // assumes no hanging nodes
        for (auto &b_vert : c_mesh->get_vertices())
        {
            if ((new_vert_2_coord - b_vert->coords()).norm() < 1E-15)
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

            c_mesh->add_i_vertex(new_vert);

            Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(1)->coords(), new_vert_2_coord);
            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            c_mesh->add_i_edge(new_edge);

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
            if (new_edges.size() != 2)
            {
                Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(1)->coords(), new_vert_2_coord);
                CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

                c_mesh->add_i_edge(new_edge);

                new_edge->set_straight();

                new_edges.push_back(new_edge);
                c_mesh->add_edge(new_edge);

                new_edge->add_vertex(bdry_edge->vertex(1));
                new_edge->add_vertex(new_vertices[1]);

                bdry_edge->vertex(1)->add_edge(new_edge);
                new_vertices[1]->add_edge(new_edge);

                bdry_edge->vertex(1)->add_vertex(new_vertices[1]);
                new_vertices[1]->add_vertex(bdry_edge->vertex(1));
            }
            assert(new_edges.size() == 2);
        }

        CurvedMesh::Edge *new_edge;

        if(make_straight && no_external)
        {
            Functional::Curve edge_param = Functional::StraightLine(new_vert_1_coord, new_vert_2_coord);
            new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
            new_edge->set_straight();
        }
        else
        {
            Functional::Curve edge_param = Functional::restriction(m_param, t0, t1);
            new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
        }

        // Functional::Curve edge_param = Functional::restriction(m_param, t0, t1);
        // CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);


        // t0 = edge_param.tmin;
        // t1 = edge_param.tmax;

        new_edges.push_back(new_edge);
        c_mesh->add_edge(new_edge);

        c_mesh->add_i_edge(new_edge);

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

        c_mesh->add_i_cell(new_cell);

        for (auto &edge : new_edges)
        {
            edge->add_cell(new_cell);
        }
        for (auto &vert : new_vertices)
        {
            vert->add_cell(new_cell);
            new_cell->add_vertex(vert);
        }

        // if ((new_edge->parameterisation().value(t0) - new_edge->vertex(0)->coords()).norm() > 1E-15)
        // {
        //     assert((new_edge->parameterisation().value(t0) - new_edge->vertex(1)->coords()).norm() < 1E-15);
        //     std::swap(t0, t1);
        // }

        // if (std::find(t_vals.begin(), t_vals.end(), t0) == t_vals.end())
        // {
        //     t_vals.push_back(t0);
        //     verts.push_back(new_edge->vertex(0));
        // }
        // if (std::find(t_vals.begin(), t_vals.end(), t1) == t_vals.end())
        // {
        //     t_vals.push_back(t1);
        //     verts.push_back(new_edge->vertex(1));
        // }

        // c_mesh->remove_edge(new_edge);
        curved_edges.push_back(new_edge);
    }

    if (no_external)
    {
        // std::cout << "hello\n";
        current_edges = c_mesh->get_edges();
        for (auto &bdry_edge : current_edges)
        {
            if (bdry_edge->n_cells() != 1)
            {
                continue;
            }

            // create two new vertices
            std::vector<CurvedMesh::Vertex *> new_vertices;
            std::vector<CurvedMesh::Edge *> new_edges;

            double t0;
            if(make_straight)
            {
                t0 = Math::atan2(bdry_edge->vertex(0)->coords()(1), bdry_edge->vertex(0)->coords()(0), 1);
            }
            else
            {
                t0 = bdry_edge->parameterisation().tmin;
            }
            Eigen::Vector2d new_vert_1_coord(std::cos(t0), std::sin(t0));

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

                c_mesh->add_i_vertex(new_vert);

                Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(0)->coords(), new_vert_1_coord);
                CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

                c_mesh->add_i_edge(new_edge);

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

            double t1;
            if(make_straight)
            {
                t1 = Math::atan2(bdry_edge->vertex(1)->coords()(1), bdry_edge->vertex(1)->coords()(0), 1);
            }
            else
            {
                t1 = bdry_edge->parameterisation().tmax;
            }
            Eigen::Vector2d new_vert_2_coord(std::cos(t1), std::sin(t1));

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

                c_mesh->add_i_vertex(new_vert);

                Functional::Curve edge_param = Functional::StraightLine(bdry_edge->vertex(1)->coords(), new_vert_2_coord);
                CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

                c_mesh->add_i_edge(new_edge);

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

            if (t0 > t1)
            {
                std::swap(t0, t1);
            }
            bool tmax_fourth_quad = (t1 > 3.0 * Math::PI / 2.0);
            bool tmin_first_quad = (t0 < Math::PI / 2.0);
            if (tmin_first_quad && tmax_fourth_quad)
            {
                t0 += 2.0 * Math::PI;
                std::swap(t0, t1);
            }

            std::function<Eigen::Vector2d(double)> circle_val_outer = [](double t) -> Eigen::Vector2d
            { return Eigen::Vector2d(std::cos(t), std::sin(t)); };
            std::function<Eigen::Vector2d(double)> circle_deriv_outer = [](double t) -> Eigen::Vector2d
            { return Eigen::Vector2d(-std::sin(t), std::cos(t)); };

            PolyMesh2D::Functional::Curve edge_param(t0, t1, circle_val_outer, circle_deriv_outer);

            CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

            // Functional::Curve edge_param = Functional::StraightLine(new_vert_1_coord, new_vert_2_coord);
            // CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
            // new_edge->set_straight();

            new_edges.push_back(new_edge);
            c_mesh->add_edge(new_edge);

            c_mesh->add_i_edge(new_edge);

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
            c_mesh->add_i_cell(new_cell);

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
    else
    {
        std::vector<CurvedMesh::Cell *> candidate_cells;

        std::vector<CurvedMesh::Vertex *> verts;
        std::vector<double> t_vals;

        for (size_t iE = 0; iE < curved_edges.size(); ++iE)
        {
            double t0 = curved_edges[iE]->parameterisation().tmin;
            double t1 = curved_edges[iE]->parameterisation().tmax;

            if ((curved_edges[iE]->parameterisation().value(t0) - curved_edges[iE]->vertex(0)->coords()).norm() > 1E-15)
            {
                assert((curved_edges[iE]->parameterisation().value(t0) - curved_edges[iE]->vertex(1)->coords()).norm() < 1E-15);
                std::swap(t0, t1);
            }

            if (std::find(t_vals.begin(), t_vals.end(), t0) == t_vals.end())
            {
                t_vals.push_back(t0);
                verts.push_back(curved_edges[iE]->vertex(0));
            }
            if (std::find(t_vals.begin(), t_vals.end(), t1) == t_vals.end())
            {
                t_vals.push_back(t1);
                verts.push_back(curved_edges[iE]->vertex(1));
            }
            for (auto &cell : curved_edges[iE]->get_cells())
            {
                candidate_cells.push_back(cell);
            }
            c_mesh->remove_edge(curved_edges[iE]);
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

        assert(verts[0] == verts[t_vals.size() - 1]);

        for (size_t i = 0; i < t_vals.size() - 1; ++i)
        {
            size_t i_plus = (i + 1) % t_vals.size();
            // Functional::Curve edge_param = Functional::restriction(m_param, t_vals[i], t_vals[i_plus]);
            CurvedMesh::Edge *new_edge;

            // Functional::Curve edge_param = Functional::StraightLine(verts[i]->coords(), verts[i + 1]->coords());
            // CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
            // new_edge->set_straight();

            if(make_straight)
            {
                Functional::Curve edge_param = Functional::StraightLine(verts[i]->coords(), verts[i + 1]->coords());
                new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
                new_edge->set_straight();
            }
            else
            {
                Functional::Curve edge_param = Functional::restriction(m_param, t_vals[i], t_vals[i_plus]);
                new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
            }

            c_mesh->add_edge(new_edge);
            c_mesh->add_i_edge(new_edge);

            new_edge->add_vertex(verts[i]);
            new_edge->add_vertex(verts[i_plus]);

            verts[i]->add_edge(new_edge);
            verts[i_plus]->add_edge(new_edge);

            verts[i]->add_vertex(verts[i_plus]);
            verts[i_plus]->add_vertex(verts[i]);

            std::vector<CurvedMesh::Cell *> joint_cells;
            for (auto &cell : candidate_cells)
            {
                std::vector<double> ts;
                for (auto &vert : cell->get_vertices())
                {
                    if (std::abs(m_level_set.value(vert->coords())) < 1E-15)
                    {
                        double t = Math::atan2(vert->coords()(1), vert->coords()(0), 1);
                        if (std::abs(t) < 1E-15)
                        {
                            if (cell->center_mass()(1) < 0.0)
                            {
                                t += 2.0 * Math::PI;
                            }
                        }
                        ts.push_back(t);
                    }
                }
                double tmin = ts[0];
                double tmax = ts[0];
                for (size_t it = 1; it < ts.size(); ++it)
                {
                    double t = ts[it];
                    if (t < tmin)
                    {
                        tmin = t;
                    }
                    if (t > tmax)
                    {
                        tmax = t;
                    }
                }
                if (((t_vals[i] > tmin - 1E-14) && (t_vals[i] < tmax + 1E-14)) && ((t_vals[i_plus] > tmin - 1E-14) && (t_vals[i_plus] < tmax + 1E-14)))
                {
                    joint_cells.push_back(cell);
                }
            }
            // assert

            // std::vector<CurvedMesh::Cell *> cells_1 = verts[i]->get_cells();
            // std::vector<CurvedMesh::Cell *> cells_2 = verts[i + 1]->get_cells();

            // std::vector<CurvedMesh::Cell *> joint_cells;

            // for (auto &cell : cells_1)
            // {
            //     if (std::find(cells_2.begin(), cells_2.end(), cell) != cells_2.end())
            //     {
            //         joint_cells.push_back(cell); // in both
            //     }
            // }

            // if (joint_cells.size() != 2)
            // {
            //     std::cout << joint_cells.size() << "\n";
            //     std::cout << joint_cells[0]->center_mass() << "\n";
            // }

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

    // for (size_t iE = 0; iE < curved_edges.size(); ++iE)
    // {
    //     std::vector<CurvedMesh::Edge *> mesh_edges = c_mesh->get_edges();
    //     if(std::find(mesh_edges.begin(), mesh_edges.end(), curved_edges[iE]) == mesh_edges.end())
    //     {
    //         continue;
    //     }
    //     double i_min = curved_edges[iE]->parameterisation().tmin;
    //     double i_max = curved_edges[iE]->parameterisation().tmax;
    //     for (size_t jE = iE + 1; jE < curved_edges.size(); ++jE)
    //     {
    //         mesh_edges = c_mesh->get_edges();
    //         if(std::find(mesh_edges.begin(), mesh_edges.end(), curved_edges[iE]) == mesh_edges.end())
    //         {
    //             break;
    //         }
    //         if(std::find(mesh_edges.begin(), mesh_edges.end(), curved_edges[jE]) == mesh_edges.end())
    //         {
    //             continue;
    //         }
    //         double j_min = curved_edges[jE]->parameterisation().tmin;
    //         double j_max = curved_edges[jE]->parameterisation().tmax;
    //         if ( ((j_min > i_min - 1E-15) && (j_max < i_max + 1E-15)) || ((i_min > j_min - 1E-15) && (i_max < j_max + 1E-15)) )
    //         // if ( ((j_min >= i_min) && (j_max <= i_max)) || ((i_min >= j_min) && (i_max <= j_max)) )
    //         {
    //             std::cout << "hello!!\n";
    //             size_t inner_edge = j_min < i_min ? iE : jE;
    //             size_t outer_edge = j_min < i_min ? jE : iE;

    //             double inner_min = j_min < i_min ? i_min : j_min;
    //             double inner_max = j_min < i_min ? i_max : j_max;
    //             double outer_min = j_min < i_min ? j_min : i_min;
    //             double outer_max = j_min < i_min ? j_max : i_max;

    //             CurvedMesh::Vertex *vert1 = curved_edges[outer_edge]->vertex(0);
    //             CurvedMesh::Vertex *vert4 = curved_edges[outer_edge]->vertex(1);
    //             if ((vert1->coords() - curved_edges[outer_edge]->parameterisation().value(outer_min)).norm() > 1E-15)
    //             {
    //                 assert((vert4->coords() - curved_edges[outer_edge]->parameterisation().value(outer_min)).norm() < 1E-15);
    //                 assert((vert1->coords() - curved_edges[outer_edge]->parameterisation().value(outer_max)).norm() < 1E-15);
    //                 std::swap(vert1, vert4);
    //             }

    //             CurvedMesh::Vertex *vert2 = curved_edges[inner_edge]->vertex(0);
    //             CurvedMesh::Vertex *vert3 = curved_edges[inner_edge]->vertex(1);
    //             if ((vert2->coords() - curved_edges[inner_edge]->parameterisation().value(inner_min)).norm() > 1E-15)
    //             {
    //                 assert((vert3->coords() - curved_edges[inner_edge]->parameterisation().value(inner_min)).norm() < 1E-15);
    //                 assert((vert2->coords() - curved_edges[inner_edge]->parameterisation().value(inner_max)).norm() < 1E-15);
    //                 std::swap(vert2, vert3);
    //             }

    //             std::vector<CurvedMesh::Edge *> edges;
    //             std::vector<CurvedMesh::Vertex *> verts;
    //             edges.push_back(curved_edges[inner_edge]);
    //             if( std::abs(outer_min - inner_min) > 1E-12 )
    //             {
    //                 Functional::Curve edge_param = Functional::restriction(m_param, outer_min, inner_min);
    //                 CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
    //                 c_mesh->add_edge(new_edge);
    //                 edges.push_back(new_edge);
    //                 verts.push_back(vert2);
    //                 vert1->add_vertex(vert2);
    //                 vert2->add_vertex(vert1);
    //                 new_edge->add_vertex(vert1);
    //                 new_edge->add_vertex(vert2);
    //                 curved_edges.push_back(new_edge);
    //             }
    //             if( std::abs(outer_max - inner_max) > 1E-12 )
    //             {
    //                 Functional::Curve edge_param = Functional::restriction(m_param, inner_max, outer_max);
    //                 CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
    //                 c_mesh->add_edge(new_edge);
    //                 edges.push_back(new_edge);
    //                 verts.push_back(vert3);
    //                 vert3->add_vertex(vert4);
    //                 vert4->add_vertex(vert3);
    //                 new_edge->add_vertex(vert3);
    //                 new_edge->add_vertex(vert4);
    //                 curved_edges.push_back(new_edge);
    //             }

    //             // std::vector<CurvedMesh::Vertex *> new_verts = curved_edges[inner_edge]->get_vertices();
    //             std::vector<CurvedMesh::Cell *> cells = curved_edges[outer_edge]->get_cells();

    //             c_mesh->remove_edge(curved_edges[outer_edge]);

    //             for (auto &cell : cells)
    //             {
    //                 for (auto &edge : edges)
    //                 {
    //                     cell->add_edge(edge);
    //                     edge->add_cell(cell);
    //                 }
    //                 cell->set_edge_orientations();
    //                 for (auto &vert : verts)
    //                 {
    //                     vert->add_cell(cell);
    //                     cell->add_vertex(vert);
    //                 }
    //             }
    //             continue;
    //         }
    //     }
    // }

    // for (size_t iE = 0; iE < curved_edges.size(); ++iE)
    // {
    //     std::vector<CurvedMesh::Edge *> mesh_edges = c_mesh->get_edges();
    //     if(std::find(mesh_edges.begin(), mesh_edges.end(), curved_edges[iE]) == mesh_edges.end())
    //     {
    //         continue;
    //     }
    //     double i_min = curved_edges[iE]->parameterisation().tmin;
    //     double i_max = curved_edges[iE]->parameterisation().tmax;
    //     for (size_t jE = iE + 1; jE < curved_edges.size(); ++jE)
    //     {
    //         mesh_edges = c_mesh->get_edges();
    //         if(std::find(mesh_edges.begin(), mesh_edges.end(), curved_edges[iE]) == mesh_edges.end())
    //         {
    //             break;
    //         }
    //         if(std::find(mesh_edges.begin(), mesh_edges.end(), curved_edges[jE]) == mesh_edges.end())
    //         {
    //             continue;
    //         }
    //         double j_min = curved_edges[jE]->parameterisation().tmin;
    //         double j_max = curved_edges[jE]->parameterisation().tmax;

    //         if ( ((i_min > j_min) && (i_min < j_max)) || ((j_min > i_min) && (j_min < i_max)) )
    //         {
    //             size_t min_edge = i_min < j_min ? iE : jE;
    //             size_t max_edge = i_min < j_min ? jE : iE;

    //             double total_min = i_min < j_min ? i_min : j_min;
    //             double other_min = i_min < j_min ? j_min : i_min;
    //             double total_max = i_min < j_min ? j_max : i_max;
    //             double other_max = i_min < j_min ? i_max : j_max;

    //             CurvedMesh::Vertex *vert1 = curved_edges[min_edge]->vertex(0);
    //             CurvedMesh::Vertex *vert3 = curved_edges[min_edge]->vertex(1);
    //             if ((vert1->coords() - curved_edges[min_edge]->parameterisation().value(total_min)).norm() > 1E-15)
    //             {
    //                 assert((vert3->coords() - curved_edges[min_edge]->parameterisation().value(total_min)).norm() < 1E-15);
    //                 std::swap(vert1, vert3);
    //             }

    //             CurvedMesh::Vertex *vert2 = curved_edges[max_edge]->vertex(0);
    //             CurvedMesh::Vertex *vert4 = curved_edges[max_edge]->vertex(1);
    //             if ((vert2->coords() - curved_edges[max_edge]->parameterisation().value(total_max)).norm() > 1E-15)
    //             {
    //                 assert((vert4->coords() - curved_edges[max_edge]->parameterisation().value(total_max)).norm() < 1E-15);
    //                 std::swap(vert2, vert4);
    //             }

    //                 vert1->add_vertex(vert2);
    //                 vert2->add_vertex(vert1);

    //                 vert2->add_vertex(vert3);
    //                 vert3->add_vertex(vert2);

    //                 vert3->add_vertex(vert4);
    //                 vert4->add_vertex(vert3);

    //             std::vector<CurvedMesh::Vertex *> verts = {vert1, vert2, vert3, vert4};

    //             Functional::Curve edge_param_1 = Functional::restriction(m_param, total_min, other_min);
    //             CurvedMesh::Edge *edge_1 = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param_1);

    //             Functional::Curve edge_param_2 = Functional::restriction(m_param, other_min, other_max);
    //             CurvedMesh::Edge *edge_2 = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param_2);

    //             Functional::Curve edge_param_3 = Functional::restriction(m_param, other_max, total_max);
    //             CurvedMesh::Edge *edge_3 = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param_3);

    //             curved_edges.push_back(edge_1);
    //             curved_edges.push_back(edge_2);
    //             curved_edges.push_back(edge_3);
    //             c_mesh->add_edge(edge_1);
    //             c_mesh->add_edge(edge_2);
    //             c_mesh->add_edge(edge_3);

    //             assert(curved_edges[min_edge]->n_cells() == 1);
    //             assert(curved_edges[max_edge]->n_cells() == 1);

    //             CurvedMesh::Cell * min_cell = curved_edges[min_edge]->cell(0);
    //             CurvedMesh::Cell * max_cell = curved_edges[max_edge]->cell(0);

    //             min_cell->add_edge(edge_1);
    //             min_cell->add_edge(edge_2);
    //             edge_1->add_cell(min_cell);
    //             edge_2->add_cell(min_cell);

    //             max_cell->add_edge(edge_2);
    //             max_cell->add_edge(edge_3);
    //             edge_2->add_cell(max_cell);
    //             edge_3->add_cell(max_cell);

    //             c_mesh->remove_edge(curved_edges[min_edge]);
    //             c_mesh->remove_edge(curved_edges[max_edge]);

    //             min_cell->add_vertex(vert2);
    //             vert2->add_cell(min_cell);

    //             max_cell->add_vertex(vert3);
    //             vert3->add_cell(max_cell);

    //             min_cell->set_edge_orientations();
    //             max_cell->set_edge_orientations();

    //         }
    //     }
    // }

    // exit(0);

    // std::cout << "\n"
    //           << verts.size() << "\n";

    // assert(verts.size() == t_vals.size());
    // assert(verts.size() > 0);

    // for (size_t i = 0; i < t_vals.size(); ++i)
    // {
    //     for (size_t j = i; j < t_vals.size(); ++j)
    //     {
    //         // if ( ((t_vals[i] < 2.0 * Math::PI) && (t_vals[j] < 2.0 * Math::PI)) || ((t_vals[i] > 2.0 * Math::PI) && (t_vals[j] > 2.0 * Math::PI)) )
    //         // {
    //             if (t_vals[i] > t_vals[j])
    //             {
    //                 std::swap(t_vals[i], t_vals[j]);
    //                 std::swap(verts[i], verts[j]);
    //             }
    //         // }
    //         // else if ( (t_vals[i] < 2.0 * Math::PI) && (t_vals[j] > 2.0 * Math::PI) )
    //         // {
    //         //     if (t_vals[i] > t_vals[j] - 2.0 * Math::PI)
    //         //     {
    //         //         std::swap(t_vals[i], t_vals[j]);
    //         //         std::swap(verts[i], verts[j]);
    //         //     }
    //         // }
    //         // else if ( (t_vals[i] > 2.0 * Math::PI) && (t_vals[j] < 2.0 * Math::PI) )
    //         // {
    //         //     if (t_vals[i] - 2.0 * Math::PI > t_vals[j])
    //         //     {
    //         //         std::swap(t_vals[i], t_vals[j]);
    //         //         std::swap(verts[i], verts[j]);
    //         //     }
    //         // }

    //     }
    // }

    // // for (size_t i = 0; i < t_vals.size(); ++i)
    // // {
    // //     std::cout << t_vals[i] << "\n";
    // // }

    // // std::cout << "\n"
    // //           << verts[0]->n_cells() << "\n";

    // // std::cout << "\n"
    // //           << verts[t_vals.size() - 1]->n_cells() << "\n";

    // // std::cout << t_vals.size() << "\n";

    // for (size_t i = 0; i < verts.size(); ++i)
    // {
    //     std::vector<CurvedMesh::Cell *> cells_i = verts[i]->get_cells();
    //     // std::vector<CurvedMesh::Cell *> cells_j = verts[j]->get_cells();
    //     for (auto &cell : cells_i)
    //     {
    //         std::vector<CurvedMesh::Vertex *> verts_on_arc;
    //         for (auto &vert : cell->get_vertices())
    //         {
    //             if (std::abs(m_level_set.value(vert->coords())) < 1E-15)
    //             {
    //                 verts_on_arc.push_back(vert);
    //             }
    //         }

    //         // double tmin = Math::atan2(verts_on_arc[0]->coords()(1), verts_on_arc[0]->coords()(0), 1);
    //         double tmin = -100.0;
    //         for(size_t iV = 0; iV < verts.size(); ++iV)
    //         {
    //             if(verts_on_arc[0] == verts[iV])
    //             {
    //                 tmin = t_vals[iV];
    //                 break;
    //             }
    //         }
    //         assert(tmin > -50.0);
    //         // double tmin = verts_on_arc[0]->parameterisation().tmin;
    //         // if()
    //         // if ((verts_on_arc[0]->parameterisation().value(t0) - new_edge->vertex(0)->coords()).norm() > 1E-15)
    //         // {
    //         //     assert((new_edge->parameterisation().value(t0) - new_edge->vertex(1)->coords()).norm() < 1E-15);
    //         //     std::swap(t0, t1);
    //         // }
    //         double tmax = tmin;

    //         for (size_t it = 1; it < verts_on_arc.size(); ++it)
    //         {
    //             // double t = std::atan2(verts_on_arc[it]->coords()(1), verts_on_arc[it]->coords()(0));
    //             double t = -100.0;
    //             for(size_t iV = 0; iV < verts.size(); ++iV)
    //             {
    //                 if(verts_on_arc[it] == verts[iV])
    //                 {
    //                     t = t_vals[iV];
    //                     break;
    //                 }
    //             }
    //             assert(t > -50.0);
    //         // double t = Math::atan2(verts_on_arc[it]->coords()(1), verts_on_arc[it]->coords()(0), 1);
    //             if (t < tmin)
    //             {
    //                 tmin = t;
    //             }
    //             if (t > tmax)
    //             {
    //                 tmax = t;
    //             }
    //         }

    //         for (size_t j = i+1; j < verts.size(); ++j)
    //         {
    //             if (std::find(verts_on_arc.begin(), verts_on_arc.end(), verts[j]) != verts_on_arc.end())
    //             {
    //                 continue;
    //             }

    //             // double t = std::atan2(verts[j]->coords()(1), verts[j]->coords()(0));
    //             double t = t_vals[j];
    //         // double t = Math::atan2(verts[j]->coords()(1), verts[j]->coords()(0), 1);

    //             if (t > tmin && t < tmax)
    //             {
    //                 verts[j]->add_cell(cell);
    //                 cell->add_vertex(verts[j]);
    //             }
    //         }
    //     }
    // }

    // for (size_t i = 0; i < (t_vals.size() - 1); ++i)
    // {
    //     Functional::Curve edge_param = Functional::restriction(m_param, t_vals[i], t_vals[i + 1]);
    //     CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);

    //     // Functional::Curve edge_param = Functional::StraightLine(verts[i]->coords(), verts[i + 1]->coords());
    //     // CurvedMesh::Edge *new_edge = new CurvedMesh::Edge(c_mesh->n_edges(), edge_param);
    //     // new_edge->set_straight();

    //     c_mesh->add_edge(new_edge);

    //     new_edge->add_vertex(verts[i]);
    //     new_edge->add_vertex(verts[i + 1]);

    //     verts[i]->add_edge(new_edge);
    //     verts[i + 1]->add_edge(new_edge);

    //     verts[i]->add_vertex(verts[i + 1]);
    //     verts[i + 1]->add_vertex(verts[i]);

    //     std::vector<CurvedMesh::Cell *> cells_1 = verts[i]->get_cells();
    //     std::vector<CurvedMesh::Cell *> cells_2 = verts[i + 1]->get_cells();

    //     std::vector<CurvedMesh::Cell *> joint_cells;

    //     for (auto &cell : cells_1)
    //     {
    //         if (std::find(cells_2.begin(), cells_2.end(), cell) != cells_2.end())
    //         {
    //             joint_cells.push_back(cell); // in both
    //         }
    //     }

    //     if (joint_cells.size() != 2)
    //     {
    //         std::cout << joint_cells.size() << "\n";
    //         std::cout << cells_1.size() << "\n";
    //         std::cout << cells_2.size() << "\n\n";
    //     }

    //     // assert(joint_cells.size() == 2); // should find exactly two cells

    //     for (auto &cell : joint_cells)
    //     {
    //         cell->add_edge(new_edge);
    //         new_edge->add_cell(cell);
    //         cell->set_edge_orientations();
    //     }
    // }
}

void MeshIntersect::make_isotropic(CurvedMesh::Mesh *c_mesh)
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

        double eps = 0.1;
        double isotropy = cell_measure / bdry_measure;
        if (isotropy < eps * diam)
        {
            CurvedMesh::Edge *longest_edge = (cell->edge(0)->n_cells() == 1 ? cell->edge(1) : cell->edge(0));
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

void MeshIntersect::make_boundary(CurvedMesh::Mesh *c_mesh)
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
