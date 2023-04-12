# CMakeList.txt : CMake project for MeshFix, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
project ("MeshFix" C CXX)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_PREFIX_PATH "D:/DEP")
set(Eigen_dir "D:/DEP/share/eigen3/cmake")
set(EIGEN3_INCLUDE_DIR "D:/DEP/include/eigen3")
find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

add_executable (MeshFix "MeshFix.cpp" "MeshFix.h" "Polyhedron.cpp"  "Polyhedron.h")
target_compile_definitions(MeshFix PRIVATE CGAL_PMP_USE_CERES_SOLVER _USE_MATH_DEFINES)
target_link_libraries(MeshFix CGAL::CGAL)
target_include_directories(MeshFix PRIVATE ${EIGEN3_INCLUDE_DIR})

# TODO: Add tests and install targets if needed.
