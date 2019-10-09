#define OM_STATIC_BUILD
#include "MIPS.h"
#include "meshio.h"

int main(int argc, char *argv[]) {
    std::string path;
    if (argc < 2) {
        std::cout << "Please specify the file name" << std::endl;
        std::cin >> path;
    }
    else path = argv[1];
    // Import the original mesh
    auto mesh = getOriginalMesh(path);
    // parameterize the mesh    
    auto param = bijectiveParameterization(mesh);
    // export the parameterization
    auto omesh = getParameterizedMesh(mesh, param);
    OpenMesh::IO::write_mesh(omesh, path + "_param.obj");
    return 0;
}