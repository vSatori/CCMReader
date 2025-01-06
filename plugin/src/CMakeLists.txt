set(classes
  vtkCCMReader)

set(CCMIO_DIR CACHE PATH "")
include_directories(${CCMIO_DIR}/include)
link_directories(${CCMIO_DIR}/lib)
set(CCMIO_LIB ccmio.lib)

vtk_module_add_module(VTK::IOCCM CLASSES ${classes})
paraview_add_server_manager_xmls(XMLS CCMReader.xml)

target_link_libraries(IOCCM PRIVATE ${CCMIO_LIB})