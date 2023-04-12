#include "MeshFix.h"
#include <atomic>
#include <cctype>
#include <list>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

#include <filesystem>
#include <utility>
#include <omp.h>
#include "happly.h"

namespace 
{
bool gVerbose = false;
}
std::pair<std::vector<Point_3>, std::vector<Triangle>> LoadOBJVF( const std::string& path )
{

    std::vector<Point_3> vertices;
    std::vector<Triangle> faces;

    std::ifstream ifs(path);
    if(!ifs.is_open())
    {
        std::cerr << "Cannot open file: " << path << std::endl;
        exit(-1);
    }
    std::string buf;
    while(std::getline(ifs, buf))
    {
        auto pos = std::find_if_not(buf.begin(), buf.end(), [](char c){ return c == ' '; });
        if(pos == buf.end())
            continue;
        if(*pos == 'v' && *std::next(pos) == ' ')
        {
            std::stringstream ss(buf);
            ss.seekg(2);
            float pos[3];
            ss >> pos[0] >> pos[1] >> pos[2];
            vertices.emplace_back(pos[0], pos[1], pos[2]);
        }
        else if (*pos == 'f')
        {
            std::stringstream ss(buf);
            ss.seekg(2);
            size_t idx[3];
            ss >> idx[0];
            ss.ignore(9999, ' ');
            ss >> idx[1];
            ss.ignore(9999, ' ');
            ss >> idx[2];
            faces.emplace_back(idx[0] - 1, idx[1] - 1, idx[2] - 1);
        }
    }

    if(gVerbose)
    {
        std::cout << "Load " << vertices.size() << " vertices," << faces.size() << " faces" << std::endl;
    }

    return {vertices, faces};
}

std::pair<std::vector<Point_3>, std::vector<Triangle>> LoadPLYVF( const std::string& path )
{
    happly::PLYData mesh(path);
    std::vector<std::array<double, 3>> ply_vertices = mesh.getVertexPositions();
    std::vector<std::vector<size_t>> ply_faces = mesh.getFaceIndices<size_t>();
    std::vector<Point_3> vertices;
    std::vector<Triangle> faces;

    for(const auto& v : ply_vertices)
    {
        vertices.emplace_back(v[0], v[1], v[2]);
    }

    for(const auto& f : ply_faces)
    {
        faces.emplace_back(f[0], f[1], f[2]);
    }

    if(gVerbose)
    {
        std::cout << "Load " << vertices.size() << " vertices," << faces.size() << " faces" << std::endl;
    }
    
    return {vertices, faces};
}

void WriteOBJVF( const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces, const std::string& path )
{
    std::ofstream ofs(path);
    for(auto& v : vertices)
    {
        ofs << "v " << v.x() << ' ' << v.y() << ' ' << v.z() << '\n';
    }

    for(auto& f : faces)
    {
        ofs << "f " << f[0] + 1 << "// " << f[1] + 1 << "// " << f[2] + 1<< "// \n";
    }

    ofs.close();
}

// void WritePLYVF( const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces, const std::string& path, int format/*0=ASCII,1=Binary,2=BinaryBigEndian*/)
// {
//     std::vector<std::array<double, 3>> ply_vertices;
//     std::vector<std::vector<size_t>> ply_faces;

//     for(int i = 0; i < vertices.size(); i++)
//     {
//         ply_vertices.emplace_back(vertices[i].x(), vertices[i].y(), vertices[i].z());
//     }

//     for(int i = 0; i < faces.size(); i++)
//     {
//         ply_faces.emplace_back( faces[i][0], faces[i][1], faces[i][2]);
//     }

//     happly::PLYData mesh;
//     mesh.addVertexPositions(ply_vertices);
//     mesh.addFaceIndices(ply_faces);

//     switch(format)
//     {
//     case 0:
//         mesh.write(path, happly::DataFormat::ASCII);
//         break;
//     case 1:
//         mesh.write(path, happly::DataFormat::Binary);
//         break;
//     case 2:
//         mesh.write(path, happly::DataFormat::BinaryBigEndian);
//         break;
//     default:
//         mesh.write(path, happly::DataFormat::ASCII);
//         break;
//     }
// }

// void WriteSTLVF( const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces, const std::string& path, int format/*0=ASCII,1=Binary*/)
// {
//     microstl::Mesh mesh;
//     for(int i = 0; i < faces.size(); i++)
//     {
//         const auto& p0 = vertices[faces[i][0]];
//         const auto& p1 = vertices[faces[i][1]];
//         const auto& p2 = vertices[faces[i][2]];
//         microstl::Facet f;
//         f.v1 = microstl::Vertex{static_cast<float>(p0.x()), static_cast<float>(p0.y()), static_cast<float>(p0.z())};
//         f.v2 = microstl::Vertex{static_cast<float>(p1.x()), static_cast<float>(p1.y()), static_cast<float>(p1.z())};
//         f.v3 = microstl::Vertex{static_cast<float>(p2.x()), static_cast<float>(p2.y()), static_cast<float>(p2.z())};
//         mesh.facets.push_back(f);
//     }
//     microstl::MeshProvider hmesh(mesh);
//     hmesh.ascii = format == 0;
//     microstl::Writer::writeStlFile(path, hmesh);
// }

std::vector<Triangle> RemoveNonManifold(const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces)
{
    std::vector<std::pair<Triangle, bool>> faceflags;
    for(auto& f : faces)
    {
        faceflags.push_back(std::make_pair(f, true));
    }

    std::unordered_map<std::pair<size_t, size_t>, Edge, PairHash, PairPred> edges;
    for(size_t i = 0; i < faceflags.size(); i++)
    {
        const auto& f = faceflags[i].first;
        auto ie0 = edges.find(std::make_pair(f[0], f[1]));
        if(ie0 == edges.end())
        {
            edges[{f[0], f[1]}] = Edge(f[0], f[1]);
            edges[{f[0], f[1]}]._faces.push_back(i);
        }
        else
        {
            edges[{f[0], f[1]}]._faces.push_back(i);
        }

        auto ie1 = edges.find({f[1], f[2]});
        if(ie1 == edges.end())
        {
            edges[{f[1], f[2]}] = Edge(f[1], f[2]);
            edges[{f[1], f[2]}]._faces.push_back(i);
        }
        else
        {
            edges[{f[1], f[2]}]._faces.push_back(i);
        }

        auto ie2 = edges.find({f[2], f[0]});
        if(ie2 == edges.end())
        {
            edges[{f[2], f[0]}] = Edge(f[2], f[0]);
            edges[{f[2], f[0]}]._faces.push_back(i);
        }
        else
        {
            edges[{f[2], f[0]}]._faces.push_back(i);
        }
    }
    
    std::vector<size_t> problematic_vertices;
    size_t nb_nm_edges = 0;
    for(auto it = edges.begin(); it != edges.end(); it++)
    {
        if(it->second._faces.size() <= 2)
        {
            continue;
        }

        problematic_vertices.push_back(it->first.first);
        problematic_vertices.push_back(it->first.first);

        for(const auto& hf : it->second._faces)
        {
            nb_nm_edges++;
            faceflags[hf].second = false;
        }
    }

    std::vector<std::vector<size_t>> vneighbors;
    vneighbors.resize(vertices.size());
    for(size_t i = 0; i < faceflags.size(); i++)
    {
        if(faceflags[i].second)
        {
            vneighbors[faceflags[i].first[0]].push_back(i);
            vneighbors[faceflags[i].first[1]].push_back(i);
            vneighbors[faceflags[i].first[2]].push_back(i);
        }
    }

    for(size_t pv : problematic_vertices)
    {
        for(size_t f : vneighbors[pv])
        {
            faceflags[f].second = false;
        }
    }

    std::atomic_int nb_nm_vertices = 0;
#pragma omp parallel for
    for(size_t iv = 0; iv < vneighbors.size(); iv++)
    {
        auto& neighbors = vneighbors[iv];
        std::list<std::pair<size_t, size_t>> sur_edges;
        size_t nb_connect_faces = neighbors.size();
        std::vector<int> sampled(nb_connect_faces, 0);
        size_t nb_cluster = 0;
        for(size_t i = 0; i < nb_connect_faces; i++)
        {
            if(sampled[i] == 1)
                continue;
            std::list<size_t> cluster;
            cluster.push_back(i);
            sampled[i] = 1;
            do
            {
                auto e0 = faceflags[neighbors[cluster.front()]].first.GetEdge(0);
                auto e1 = faceflags[neighbors[cluster.front()]].first.GetEdge(1);
                auto e2 = faceflags[neighbors[cluster.front()]].first.GetEdge(2);

                for(size_t j = 0; j < nb_connect_faces; j++)
                {
                    if(j != cluster.front() && sampled[j] != 1)
                    {
                        auto e3 = faceflags[neighbors[j]].first.GetEdge(0);
                        auto e4 = faceflags[neighbors[j]].first.GetEdge(1);
                        auto e5 = faceflags[neighbors[j]].first.GetEdge(2);

                        if(PairPred()(e0, e3) || PairPred()(e0, e4) || PairPred()(e0, e5) ||
                        PairPred()(e1, e3) || PairPred()(e1, e4) || PairPred()(e1, e5) ||
                        PairPred()(e2, e3) || PairPred()(e2, e4) || PairPred()(e2, e5))
                        {
                            cluster.push_back(j);
                            sampled[j] = 1;
                        }
                    }
                }
                cluster.pop_front();
            } while(!cluster.empty());
            nb_cluster++;
        }

        if(nb_cluster > 1)
        {
            nb_nm_vertices++;
            for(size_t hf : neighbors)
            {
                faceflags[hf].second = false;
            }
        }
    }

    std::vector<Triangle> result_faces;
    for(const auto& [face, flag] : faceflags)
    {
        if(flag)
        {
            result_faces.push_back(face);
        }
    }

    if(gVerbose)
    {
        std::cout << "Find " << nb_nm_edges << " non-manifold edges and " << nb_nm_vertices << " non-manifold vertices." << std::endl;
        std::cout << "After remove non-manifold: " << result_faces.size() << " faces." << std::endl;
    }
    return result_faces;
}

bool IsSmallHole( hHalfedge hh, Polyhedron& mesh, int max_num_hole_edges, float max_hole_diam)
{
    int num_hole_edges = 0;
    CGAL::Bbox_3 hole_bbox;
    for (hHalfedge hc : CGAL::halfedges_around_face(hh, mesh))
    {
        const Point_3& p = hc->vertex()->point();
        hole_bbox += p.bbox();
        ++num_hole_edges;
        // Exit early, to avoid unnecessary traversal of large holes
        if (num_hole_edges > max_num_hole_edges) return false;
        if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
        if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
        if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
    }
    return true;
}

int main(int argc, char* argv[])
{
    auto print_help_msg = []()
    {
        std::cout << "usage:\n"
        "\t-i filename \tPath to input mesh. OBJ format only.\n"
        "\t-o filename \tFile name of output mesh. OBJ format only.\n"
        "\t-k \tKeep largest connected component (default=off)\n"
        "\t-s \tFix self intersection\n"
        "\t-f max_hole_edges max_hole_diam\t Do not fill big holes that satisfiy (edge number > max_hole_edges) or (bounding box size > max_hole_diam)\n"
        "\t-r refine holes.\n"
        "\t-v \tPrint debug messages (default=off)" << std::endl;
    };
    if(argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
    {
        print_help_msg();
        return -1;
    }
    std::string path;
    std::string output_path;
    bool keep_largest_connected_component = false;
    bool fix_self_intersection = false;
    bool filter_small_holes = false;
    int max_hole_edges = std::numeric_limits<int>::max();
    float max_hole_diam = std::numeric_limits<float>::max();
    bool refine = false;

    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-i") == 0)
        {
            path = std::string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-o") == 0)
        {
            output_path = std::string(argv[i + 1]);
        }
        else if (strcmp(argv[i], "-v") == 0)
        {
            gVerbose = true;
        }
        else if (strcmp(argv[i], "-k") == 0)
        {
            keep_largest_connected_component = true;
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            fix_self_intersection = true;
        }
        else if (strcmp(argv[i], "-f") == 0)
        {
            filter_small_holes = true;
            max_hole_edges = std::atoi(argv[i+1]);
            max_hole_diam = std::atof(argv[i+2]);
        }
    }
    if(path.empty() || output_path.empty())
    {
        print_help_msg();
        return -1;
    }
    
    std::string postfix = path.substr(path.rfind('.'));
    std::pair<std::vector<Point_3>, std::vector<Triangle>> pair;
    if(postfix == ".obj")
        pair = LoadOBJVF(path);
    else if (postfix == ".ply")
        pair = LoadPLYVF(path);
    else
    {
        print_help_msg();
        return -1;
    }
    std::vector<Point_3> vertices = std::move(pair.first);
    std::vector<Triangle> faces = std::move(pair.second);

    auto new_faces = RemoveNonManifold(vertices, faces);

    std::vector<int> indices;
    for(const auto& f : new_faces )
    {
        indices.push_back(f[0]);
        indices.push_back(f[1]);
        indices.push_back(f[2]);
    }

    Polyhedron m( vertices, indices );
    
    CGAL::Polygon_mesh_processing::remove_isolated_vertices(m);

    if(fix_self_intersection)
    {
        std::vector<std::pair<hFacet, hFacet>> intersect_faces;
        CGAL::Polygon_mesh_processing::self_intersections<CGAL::Parallel_if_available_tag>(m, std::back_inserter(intersect_faces));
        std::unordered_set<hFacet> face_to_remove;
        for(auto [f1, f2] : intersect_faces)
        {
            face_to_remove.insert(f1);
            face_to_remove.insert(f2);
        }
        for(auto& hf : face_to_remove)
        {
            m.erase_facet(hf->halfedge());
        }

        auto [vertices1, faces1] = m.ToVerticesFaces();
        std::vector<Triangle> triangles1;
        for(int i = 0; i < faces1.size() / 3; i++)
        {
            triangles1.emplace_back(faces1[i * 3 + 0], faces1[i * 3 + 1], faces1[i * 3 + 2]);
        }
        auto newfaces1 = RemoveNonManifold(vertices1, triangles1);
        std::vector<int> indices1;
        for(const auto& f : newfaces1 )
        {
            indices1.push_back(f[0]);
            indices1.push_back(f[1]);
            indices1.push_back(f[2]);
        }
        
        m = Polyhedron(vertices1, indices1);
    }

    if(keep_largest_connected_component)
    {
        CGAL::Polygon_mesh_processing::keep_largest_connected_components(m, 1);
    }

    std::vector<hHalfedge> border_edges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m, std::back_inserter(border_edges));
    for(hHalfedge hh : border_edges)
    {
        if(filter_small_holes)
        {
            if(IsSmallHole(hh, m, max_hole_edges, max_hole_diam))
            {
                if(refine)
                {
                    std::vector<hVertex> patch_vertices;
                    std::vector<hFacet> patch_faces;
                    CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(m, hh, std::back_inserter(patch_faces), std::back_inserter(patch_vertices));
                }
                else
                {
                    std::vector<hFacet> patch_faces;
                    CGAL::Polygon_mesh_processing::triangulate_hole(m, hh, std::back_inserter(patch_faces));
                }
            }
        }
        else
        {
            if(refine)
            {
                std::vector<hVertex> patch_vertices;
                std::vector<hFacet> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(m, hh, std::back_inserter(patch_faces), std::back_inserter(patch_vertices));
            }
            else
            {
                std::vector<hFacet> patch_faces;
                CGAL::Polygon_mesh_processing::triangulate_hole(m, hh, std::back_inserter(patch_faces));
            }
        }
    }
    std::string out_postfix = output_path.substr(output_path.rfind('.'));
    if(out_postfix == ".ply")
        CGAL::IO::write_PLY(output_path, m);
    else if (out_postfix == ".stl")
        CGAL::IO::write_STL(output_path, m, CGAL::parameters::use_binary_mode(true));
    else
        CGAL::IO::write_OBJ(output_path, m, CGAL::parameters::use_binary_mode(true));

    if(gVerbose)
    {
        std::cout << "Output " << m.size_of_vertices() << " vertices," << m.size_of_facets() << " faces" << std::endl;
    }
    return 0;
}