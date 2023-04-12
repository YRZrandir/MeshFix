//
// Created by yrz on 7/7/22.
//

#ifndef ELASTICITY_POLYHEDRON_H
#define ELASTICITY_POLYHEDRON_H
#include <memory>
#include <utility>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

using KernelEpick = CGAL::Exact_predicates_inexact_constructions_kernel;
using CPolyhedron = CGAL::Polyhedron_3<KernelEpick>;
using hHalfedge = CPolyhedron::Halfedge_handle;
using hVertex = CPolyhedron::Vertex_handle;
using hFacet = CPolyhedron::Facet_handle;
using Halfedge = CPolyhedron::Halfedge;
using CVertex = CPolyhedron::Vertex;
using Facet = CPolyhedron::Facet;
using iHalfedge = CPolyhedron::Halfedge_iterator;
using iVertex = CPolyhedron::Vertex_iterator;
using iFacet = CPolyhedron::Facet_iterator;
using Point_3 = CPolyhedron::Point_3;
using Vec3 = CPolyhedron::Traits::Vector_3;

class Polyhedron : public CPolyhedron
{
public:
    Polyhedron( const std::vector<Point_3>& vertices, const std::vector<int>& indices );
    Polyhedron( const Polyhedron& mesh );

    std::pair<std::vector<Point_3>, std::vector<int>> ToVerticesFaces() const;

    void WriteOFF( const std::string& path ) const;
    void WriteOBJ( const std::string& path ) const;
};

template <typename HDS>
class PolyhedronObjBulider : public CGAL::Modifier_base<HDS>
{
public:
    PolyhedronObjBulider( const std::vector<Point_3>& vertices, const std::vector<int>& indices )
        : _vertices(vertices), _indices(indices) {}
    virtual void operator()( HDS& hds ) override;

protected:
    const std::vector<Point_3>& _vertices;
    const std::vector<int>& _indices;
};

template<typename HDS>
inline void PolyhedronObjBulider<HDS>::operator()( HDS& hds )
{
    CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true );
    builder.begin_surface( _vertices.size(), _indices.size() / 3 );
    for (size_t i = 0, size = _vertices.size(); i < size; i += 1)
    {
        builder.add_vertex( _vertices[i] );
    }
    for (int f = 0, size = _indices.size() / 3; f < size; ++f)
    {
        builder.begin_facet();
        builder.add_vertex_to_facet( _indices[f * 3 + 0]);
        builder.add_vertex_to_facet( _indices[f * 3 + 1] );
        builder.add_vertex_to_facet( _indices[f * 3 + 2] );
        builder.end_facet();
    }
    builder.end_surface();
}

#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME Polyhedron
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CPolyhedron
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>
#undef CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME
#undef CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME
#endif //ELASTICITY_POLYHEDRON_H
