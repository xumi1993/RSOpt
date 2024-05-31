
file(GLOB SURFKER_SOURCES "*.f90")
include_directories(
  ${PROJECT_SOURCE_DIR}/include
)
add_library(shared_obj OBJECT ${SURFKER_SOURCES})