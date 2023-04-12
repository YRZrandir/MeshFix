//
// Created by yrz on 7/7/22.
//

#include <cmath>
#include <sstream>
#include "Polyhedron.h"
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/linear_least_squares_fitting_3.h>



Polyhedron::Polyhedron(  const std::vector<Point_3>& vertices, const std::vector<int>& indices  )
{
    PolyhedronObjBulider<CPolyhedron::HalfedgeDS> builder(vertices, indices);
    delegate(builder);
}

Polyhedron::Polyhedron( const Polyhedron& mesh )
    :CPolyhedron( mesh )
{
}

std::pair<std::vector<Point_3>, std::vector<int>> Polyhedron::ToVerticesFaces() const
{
    std::vector<Point_3> vertices;
    std::vector<int> indices;

    std::unordered_map<decltype(vertices_begin()), size_t> vertex_id_map;
    size_t v_count = 0;
    for(auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        vertex_id_map[hv] = v_count++;
        vertices.push_back(hv->point());
    }

    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        size_t v0 = vertex_id_map[hf->halfedge()->vertex()];
        size_t v1 = vertex_id_map[hf->halfedge()->next()->vertex()];
        size_t v2 = vertex_id_map[hf->halfedge()->prev()->vertex()];
        indices.push_back(v0);
        indices.push_back(v1);
        indices.push_back(v2);
    }

    return {vertices, indices};
}

void Polyhedron::WriteOFF( const std::string& path ) const
{
    std::ofstream ofs( path );
    CGAL::IO::write_OFF( ofs, *this );
}

void Polyhedron::WriteOBJ( const std::string& path ) const
{
    std::stringstream ss;

    std::unordered_map<decltype(vertices_begin()), size_t> vertex_id_map;
    size_t v_count = 1;
    for(auto hv = vertices_begin(); hv != vertices_end(); hv++)
    {
        vertex_id_map[hv] = v_count++;
        ss << "v " << hv->point().x() << ' ' << hv->point().y() << ' ' << hv->point().z() << '\n';
    }

    std::unordered_map<decltype(halfedges_begin()), size_t> hedge_id_map;
    size_t h_count = 1;
    for(auto hh = halfedges_begin(); hh != halfedges_end(); hh++)
    {
        hedge_id_map[hh] = h_count++;
    }

    for(auto hf = facets_begin(); hf != facets_end(); hf++)
    {
        size_t e0 = hedge_id_map[hf->halfedge()];
        size_t e1 = hedge_id_map[hf->halfedge()->next()];
        size_t e2 = hedge_id_map[hf->halfedge()->prev()];
        size_t v0 = vertex_id_map[hf->halfedge()->vertex()];
        size_t v1 = vertex_id_map[hf->halfedge()->next()->vertex()];
        size_t v2 = vertex_id_map[hf->halfedge()->prev()->vertex()];
        ss << "f " << v0 << "//" << e0 << ' ' << v1 << "//" << e1 << ' ' << v2 << "//" << e2 << '\n';
    }

    std::ofstream ofs(path);
    ofs << ss.rdbuf();
    ofs.close();
}
