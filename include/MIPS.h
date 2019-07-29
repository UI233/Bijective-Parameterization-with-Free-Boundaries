#ifndef MIPS_H
#define MIPS_H
#include "OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh"
#include <vector>
using Mesh = OpenMesh::PolyMesh_ArrayKernelT<>;

/*
    Parameterize the mesh using MIPS
    The optimization strategy is gradient descent
    @param: Input Mesh and error bound for paratimization
    @return: the parametered point order by the order by the points of input mesh
*/
std::vector<OpenMesh::Vec2f> paratimization(Mesh&, double error = 0.0001);
#endif