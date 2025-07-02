cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(GridViewer)

list(APPEND CMAKE_PREFIX_PATH "/home/paboyle/Visualisation/install/")

find_package(VTK COMPONENTS 
  CommonColor
  CommonCore
  FiltersCore
  FiltersModeling
  IOImage
  IOFFMPEG
  InteractionStyle
  InteractionWidgets
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "GridViewer: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(FieldDensityAnimate MACOSX_BUNDLE FieldDensityAnimate.cxx )
  target_link_libraries(FieldDensityAnimate PRIVATE ${VTK_LIBRARIES}
)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS FieldDensityAnimate
  MODULES ${VTK_LIBRARIES}
)
