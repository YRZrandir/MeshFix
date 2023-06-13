# MeshFix
Fix non-manifold and self-intersection of triangular meshes.

+ The program deletes non-manifold & self-intersection parts of the mesh, and fills the holes.
+ Support format: Any format that Assimp supports.
+ Dependencies: CGAL, Assimp

# Build
Linux:
```
mkdir build
cd build
cmake ..
make
```
Windows: Use cmake-gui and it should generate a Visual Studio Project.
# Usage

+ **-i filename**  Path to input mesh.
+ **-o filename**  File name of output mesh.
+ **-k threshold** Delete connected components smaller than threshold (default=off)
+ **-s**  When turned on, the program will check self intersection.
+ **-f max_hole_edges max_hole_diam** Do not fill holes that satisfiy (edge number > max_hole_edges) OR (AABB size > max_hole_diam)
+ **-r** Refine after filling holes.
+ **-m max_retry** The program repeatedly trys to fix the mesh until no non-manifold or self-intersection is found, this is the max retry time. (default=10)
+ **-v** Print debug messages