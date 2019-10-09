#ifndef MESHIO_H
#define MESHIO_H
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#include <string>
#include <iostream>
#include <fstream>
using Mesh = OpenMesh::PolyMesh_ArrayKernelT<>;

// Read the mesh from specified file and check the constrain of mesh
// @param: the path to the file for the mesh to be parameterized
// @return: the triangular mesh obtained
// @throw: std::logic_error if the mesh has no boundary
Mesh getOriginalMesh(const std::string& path);

// Get a mesh which specifies the parameterization from original mesh and parameters for each vertices
// @param: the original mesh and the parameters for the mesh
// @return: a planar mesh specifying the parameterization
Mesh getParameterizedMesh(const Mesh& original_mesh, const std::vector<OpenMesh::Vec2f> &param);
#endif 
