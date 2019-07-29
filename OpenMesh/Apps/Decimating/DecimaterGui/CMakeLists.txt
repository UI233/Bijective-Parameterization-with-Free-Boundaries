include (ACGCommon)

include_directories (
  ../../../..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# source code directories
set (directories 
  ../../QtViewer
  ../
)

# collect all header and source files
acg_append_files (headers "*.hh" ${directories})
list (APPEND headers "../../QtViewer/MeshViewerWidgetT_impl.hh")


set (sources
  ../../QtViewer/QGLViewerWidget.cc
  ../DecimaterViewerWidget.cc
  ../decimaterviewer.cc
)

if (WIN32)
  acg_add_executable (DecimaterGui WIN32 ${sources} ${headers})
else ()
  acg_add_executable (DecimaterGui ${sources} ${headers} ) 
endif ()

target_link_libraries (DecimaterGui
  OpenMeshCore
  OpenMeshTools
  Qt5::OpenGL
  ${OPENGL_LIBRARIES}
)

