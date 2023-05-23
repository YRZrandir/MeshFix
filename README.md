# MeshFix
Fix non-manifold and self-intersection of triangular meshes.

+ The program will delete non-manifold & self-intersection part of the mesh, and fill the holes.
+ Support format: Any format that Assimp supports.
+ Dependencies: CGAL, Assimp

# Build
On Linux:
```
mkdir build
cd build
cmake ..
make
```
On Windows: Use cmake-gui and it should generate a Visual Studio Project
# Usage

+ **-i filename**  Path to input mesh.
+ **-o filename**  File name of output mesh.
+ **-k threshold** Delete connected components smaller than threshold (default=off)
+ **-s**  Fix self intersection
+ **-f max_hole_edges max_hole_diam** Do not fill holes that satisfiy (edge number > max_hole_edges) OR (AABB size > max_hole_diam)
+ **-r** refine after filling holes.
+ **-m max_retry** The program will repeatedly try to fix the mesh, this is the max retry time. (default=10)
+ **-v** Print debug messages