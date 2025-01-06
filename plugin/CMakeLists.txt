paraview_add_plugin(CCMReader
  VERSION "1.0"
  MODULES VTK::IOCCM
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/vtk.module")