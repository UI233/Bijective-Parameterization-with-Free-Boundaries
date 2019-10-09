#ifndef MIPS_H
#define MIPS_H
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#include <vector>
using Mesh = OpenMesh::PolyMesh_ArrayKernelT<>;

// Parameterize the mesh using Bijective Paramerization
// The optimization strategy is LBFGS
// @param: Input Mesh and error bound for paratimization
// @return: the parametered point order by the order by the points of input mesh
std::vector<OpenMesh::Vec2f> bijectiveParameterization(Mesh&, double error = 0.0001);
#endif