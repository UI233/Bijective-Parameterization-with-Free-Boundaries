#define OM_STATIC_BUILD
#include "MIPS.h"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
    std::string path;
    if (argc < 2) {
        std::cout << "Please specify the file name" << std::endl;
        std::cin >> path;
    }
    else path = argv[1];
    // Import mesh
    Mesh mesh;
    OpenMesh::IO::read_mesh(mesh, argv[1]);
    mesh.triangulate();
    // for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
    //     point = mesh.point(*it);
    //     std::cout << point[0] << " " << point[1] << " " << point[2] << std::endl;
    // }

    auto param = paratimization(mesh);
  
    Mesh omesh;
    for (auto v: param) {
        OpenMesh::Vec3f vert(v[0], v[1], 0.0f);
        omesh.add_vertex(vert);
    }
    for (auto face: mesh.all_faces()) {
        int cnt = 0;
        OpenMesh::VertexHandle handles[10];
        for (auto v: mesh.fv_range(face)) 
            handles[cnt++] = v;
        omesh.add_face(handles, cnt);
    }

    OpenMesh::IO::write_mesh(omesh, path + "_param.obj");

    return 0;
}