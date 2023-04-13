# CMakeList.txt : CMake project for MeshFix, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project ("MeshFix" C CXX)

set(CMAKE_BUILD_TYPE Release)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

add_executable (MeshFix "MeshFix.cpp" "MeshFix.h" "Polyhedron.cpp"  "Polyhedron.h")
target_compile_definitions(MeshFix PRIVATE _USE_MATH_DEFINES)
target_link_libraries(MeshFix CGAL::CGAL OpenMP::OpenMP_CXX)
target_include_directories(MeshFix PRIVATE "./third-party/include")
