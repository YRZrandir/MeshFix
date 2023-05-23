#include "MeshFix.h"
#include <atomic>
#include <cctype>
#include <list>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <assimp/Importer.hpp>
#include <assimp/Exporter.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <omp.h>

namespace 
{
bool gVerbose = false;
}
std::pair<std::vector<Point_3>, std::vector<Triangle>> LoadVFAssimp( const std::string& path )
{
    Assimp::Importer importer;
    importer.SetPropertyInteger(AI_CONFIG_PP_RVC_FLAGS, aiComponent_COLORS | aiComponent_NORMALS | aiComponent_TEXCOORDS );
    const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_JoinIdenticalVertices | aiProcess_RemoveComponent);
    if(!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode || scene->mNumMeshes < 1)
    {
        std::cout << "Input invalid." << std::endl;
        exit(-1);
    }
    const aiMesh* mesh = scene->mMeshes[0];

    auto to_cgal = [](const aiVector3D& v) { return Point_3{v.x, v.y, v.z};};

    std::vector<Point_3> vertices;
    for(int i = 0; i < mesh->mNumVertices; i++)
    {
        vertices.push_back(to_cgal(mesh->mVertices[i]));
    }

    std::vector<Triangle> faces;
    for(int i = 0; i < mesh->mNumFaces; i++)
    {
        const auto& f = mesh->mFaces[i];
        faces.emplace_back(f.mIndices[0], f.mIndices[1], f.mIndices[2]);
    }
    if(gVerbose)
    {
        std::cout << "Loading " << vertices.size() << " vertices, " << faces.size() << " faces. " << std::endl;
    }
    return {vertices, faces};
}

void WriteVFAssimp( const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces, const std::string& path)
{
    Assimp::Exporter exporter;
    auto scene = std::make_unique<aiScene>();

    scene->mRootNode = new aiNode();
    scene->mRootNode->mNumMeshes = 1;
    scene->mRootNode->mMeshes = new unsigned[1];
    scene->mRootNode->mMeshes[0] = 0;

    scene->mNumMaterials = 1;
    scene->mMaterials = new aiMaterial*[]{ new aiMaterial() };
    scene->mMetaData = new aiMetadata();

    scene->mNumMeshes = 1;
    scene->mMeshes = new aiMesh*[1];
    scene->mMeshes[0] = new aiMesh();

    aiMesh* m = scene->mMeshes[0];
    m->mNumFaces = faces.size();
    m->mNumVertices = vertices.size();
    m->mVertices = new aiVector3D[m->mNumVertices];
    m->mFaces = new aiFace[m->mNumFaces];
    m->mPrimitiveTypes = aiPrimitiveType_TRIANGLE;

    for(int i = 0; i < m->mNumVertices; i++)
    {
        m->mVertices[i] = aiVector3D(vertices[i].x(), vertices[i].y(), vertices[i].z());
        
    }

    for(int i = 0; i < m->mNumFaces; i++)
    {
        m->mFaces[i].mNumIndices = 3;
        m->mFaces[i].mIndices = new unsigned int[3];
        m->mFaces[i].mIndices[0] = faces[i][0];
        m->mFaces[i].mIndices[1] = faces[i][1];
        m->mFaces[i].mIndices[2] = faces[i][2];
    }

    std::string postfix = path.substr(path.rfind('.') + 1);
    
    if(postfix == std::string("ply"))
    {
        postfix = "plyb";
    }
    else if (postfix == std::string("stl"))
    {
        postfix = "stlb";
    }
    //Assimp::ExportProperties prop;
    exporter.Export(scene.get(), postfix, path);
}

void WriteCgalPolyAssimp( const Polyhedron& m, const std::string& path )
{
    auto [vertices, faces] = PolyhedronToVF( m );
    WriteVFAssimp(vertices, faces, path);
}

std::pair<std::vector<Point_3>, std::vector<Triangle>> PolyhedronToVF( const Polyhedron& m )
{
    std::vector<Point_3> vertices;

    std::unordered_map<hVertex, unsigned> idmap;
    int count = 0;
    for(auto hv : CGAL::vertices(m))
    {
        idmap[hv] = count;
        vertices.push_back(hv->point());
        count++; 
    }

    std::vector<Triangle> triangles;
    for(auto hf : CGAL::faces(m))
    {
        unsigned i0 = idmap[hf->halfedge()->vertex()];
        unsigned i1 = idmap[hf->halfedge()->next()->vertex()];
        unsigned i2 = idmap[hf->halfedge()->prev()->vertex()];
        triangles.emplace_back(i0, i1, i2);
    }

    return {vertices, triangles};
}

std::vector<Triangle> RemoveNonManifold(const std::vector<Point_3>& vertices, const std::vector<Triangle>& faces, int* nb_removed_face )
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
    for(int iv = 0; iv < vneighbors.size(); iv++)
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
    *nb_removed_face = faces.size() - result_faces.size();
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

void Run(
    std::string path,
    std::string output_path, 
    bool keep_largest_connected_component,
    int large_cc_threshold,
    bool fix_self_intersection,
    bool filter_small_holes,
    int max_hole_edges,
    float max_hole_diam,
    bool refine,
    int max_retry)
{
    auto [vertices, faces] = LoadVFAssimp(path);

    int nb_removed_faces = 0;
    int cnt = 0;
    do
    {
        faces = RemoveNonManifold(vertices, faces, &nb_removed_faces);
        if(cnt++ > max_retry)
            break;
    } while (nb_removed_faces != 0);

    std::vector<int> indices;
    for(const auto& f : faces )
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
        nb_removed_faces = 0;
        int cnt = 0;
        do
        {
            triangles1 = RemoveNonManifold(vertices1, triangles1, &nb_removed_faces);
            if(cnt++ > max_retry)
                break;
        } while(nb_removed_faces != 0);

        std::vector<int> indices1;
        for(const auto& f : triangles1 )
        {
            indices1.push_back(f[0]);
            indices1.push_back(f[1]);
            indices1.push_back(f[2]);
        }
        
        m = Polyhedron(vertices1, indices1);
    }

    if(keep_largest_connected_component)
    {
        int num = CGAL::Polygon_mesh_processing::keep_large_connected_components(m, large_cc_threshold);
        if(gVerbose)
        {
            std::cout << "Remove " << num << " small connected components." << std::endl;
        }
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

    auto [cgal_vertices, cgal_faces] = PolyhedronToVF(m);
    WriteVFAssimp(cgal_vertices, cgal_faces, output_path);

    if(gVerbose)
    {
        std::cout << "Output " << m.size_of_vertices() << " vertices," << m.size_of_facets() << " faces" << std::endl;
    }
}

int main(int argc, char* argv[])
{
    auto print_help_msg = []()
    {
        std::cout << "usage:\n"
        "\t-i filename \tPath to input mesh.\n"
        "\t-o filename \tFile name of output mesh.\n"
        "\t-k threshold\tDelete connected components smaller than threshold (default=off)\n"
        "\t-s \tFix self intersection\n"
        "\t-f max_hole_edges max_hole_diam\t Do not fill holes that satisfiy (edge number > max_hole_edges) OR (AABB size > max_hole_diam)\n"
        "\t-r refine after filling holes.\n"
        "\t-m max_retry The program will repeatedly try to fix the mesh, this is the max retry time. (default=10)"
        "\t-v \tPrint debug messages" << std::endl;
    };
    if(argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
    {
        print_help_msg();
        return -1;
    }
    std::string path;
    std::string output_path;
    bool keep_largest_connected_component = false;
    int large_cc_threshold = 100;
    bool fix_self_intersection = false;
    bool filter_small_holes = false;
    int max_hole_edges = std::numeric_limits<int>::max();
    float max_hole_diam = std::numeric_limits<float>::max();
    bool refine = false;
    int max_retry = 10;

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
            if(i < argc - 1 && std::atoi(argv[i+1]) != 0)
            {
                large_cc_threshold = std::atoi(argv[i + 1]);
            }
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
        else if (strcmp(argv[i], "-r") == 0)
        {
            refine = true;
        }
        else if (strcmp(argv[i], "-m") == 0)
        {
            max_retry = std::atoi(argv[i+1]);
        }
    }
    if(path.empty() || output_path.empty())
    {
        print_help_msg();
        return -1;
    }
    
    Run( 
        path,
        output_path,
        keep_largest_connected_component,
        large_cc_threshold,
        fix_self_intersection, 
        filter_small_holes, 
        max_hole_edges, 
        max_hole_diam, 
        refine,
        max_retry
    );
    
    return 0;
}