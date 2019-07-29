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
  ../SubdivideWidget.cc
  ../qtsubdivider.cc
)

if (WIN32)
  acg_add_executable (SubdividerGui WIN32 ${sources} ${headers})
else ()
  acg_add_executable (SubdividerGui ${sources} ${headers})
endif ()

target_link_libraries (SubdividerGui
  OpenMeshCore
  OpenMeshTools
  Qt5::OpenGL
  ${OPENGL_LIBRARIES}
)

