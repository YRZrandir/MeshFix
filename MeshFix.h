#ifndef MESH_FIX
#define MESH_FIX
#include <functional>
#include <string>
#include <utility>
#include <tuple>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include "Polyhedron.h"

struct Triangle
{
public:
    Triangle( unsigned i0, unsigned i1, unsigned i2)
    {
        _id[0] = i0;
        _id[1] = i1;
        _id[2] = i2;
    }
    size_t& operator[](size_t i) { return _id[i]; }
    const size_t& operator[](size_t i) const { return _id[i]; }
    std::pair<size_t, size_t> GetEdge(size_t i)
    {
        switch (i)
        {
        case 0:
            return std::make_pair(_id[0], _id[1]);
            break;
        case 1:
            return std::make_pair(_id[1], _id[2]);
            break;
        case 2:
            return std::make_pair(_id[2], _id[0]);
            break;
        }
        return std::make_pair<size_t, size_t>(0, 0);
    }
protected:
    size_t _id[3]{0, 0, 0};
};

struct Edge
{
public:
    Edge() = default;

    Edge(size_t i0, size_t i1)
    {
        _i0 = i0;
        _i1 = i1;
    }

    size_t _i0;
    size_t _i1;
    std::vector<size_t> _faces;
};

struct PairHash
{
    size_t operator()(const std::pair<size_t, size_t>& p) const
    {
        if(p.first > p.second)
        {
            return boost::hash<std::pair<size_t, size_t>>()(p);
        }
        else
        {
            return boost::hash<std::pair<size_t, size_t>>()({p.second, p.first});
        }
    }
};

struct PairPred
{
    bool operator()(std::pair<size_t, size_t> e0, std::pair<size_t, size_t> e1) const
    {
        return e0.first == e1.first && e0.second == e1.second || e0.first == e1.second && e0.second == e1.first;
    }
};

using GridPos = std::tuple<int, int, int>;

struct GridHash
{
    size_t operator()(const GridPos& p) const
    {
        return std::get<0>(p) * 73856093 ^ std::get<1>(p) * 19349663 ^ std::get<2>(p) * 83492791;
    }
};

struct GridPred
{
    bool operator()(const GridPos& p1, const GridPos& p2) const
    {
        return std::get<0>(p1) == std::get<0>(p2) && std::get<1>(p1) == std::get<1>(p2) && std::get<2>(p1) == std::get<2>(p2);
    }
};

class SpatialHashing
{
public:
    SpatialHashing( float dx );
    void Insert(const Triangle& t, const std::vector<Point_3>& points);
    bool CheckIntersection(const Triangle& t, const std::vector<Point_3>& points) const;
protected:
    int GridCoord( float p ) const;
protected:
    float _dx;
    std::unordered_map<GridPos, std::vector<Triangle>, GridHash, GridPred> _table;
};

std::pair<std::vector<Point_3>, std::vector<Triangle>> LoadOBJVF( const std::string& path );
std::vector<Triangle> RemoveNonManifold(const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces);
bool IsSmallHole( hHalfedge hh, Polyhedron& mesh, int max_num_hole_edges, float max_hole_diam);

#endif