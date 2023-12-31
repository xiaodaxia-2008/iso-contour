# iso-contour
Extracting ISO contour of a mesh vertex scalar field.

## Problem
Supposing there is a mesh, at each vertex there is a scalar property, we want to find the polylines where all the scalar property is constant. The polylines are a ISO contour of  this mesh. For example, if the scalar property equals to coordinate z, then the iso contour are the points whose z valves all equals to a constant value.

## Dependecy
The library uses [CGAL](https://github.com/CGAL/cgal) Halfedge structure.

## Example
The usage is simple and easy, please refer to [geodesic_example](./examples/geodesic_contour.cpp)
![sdf](./doc/sdf.png)

## Building

You can use `vcpkg` to install `cgal`， and then build and install with normal cmake procedure:

```shell
cd iso-contour
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
cmake --install . --prefix /path/prefix/to/install
```