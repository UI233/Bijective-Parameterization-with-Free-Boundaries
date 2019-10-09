#include "meshio.h"
#include <exception>

Mesh getOriginalMesh(const std::string& path) {
    // read and triangulate the mesh
    Mesh mesh;
    OpenMesh::IO::read_mesh(mesh, path);
    mesh.triangulate();
    
    // check the existance of boudary
    for (auto vitr = mesh.vertices_begin(); vitr != mesh.vertices_end(); vitr++) {
        if (mesh.is_boundary(*vitr))
            return mesh;
    }

    throw std::logic_error("The mesh has no boudary");
}

Mesh getParameterizedMesh(const Mesh& original_mesh, const std::vector<OpenMesh::Vec2f> &param) {
    Mesh omesh;
    for (auto v: param) {
        OpenMesh::Vec3f vert(v[0], v[1], 0.0f);
        omesh.add_vertex(vert);
    }
    for (auto face: original_mesh.all_faces()) {
        int cnt = 0;
        OpenMesh::VertexHandle handles[10];
        for (auto v: original_mesh.fv_range(face)) 
            handles[cnt++] = v;
        omesh.add_face(handles, cnt);
    }
    return omesh;
}