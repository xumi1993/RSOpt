# fortran-yaml-cpp
add_subdirectory("${SOURCE_DIR}/fortran-yaml")

# Shared library
add_subdirectory("${SOURCE_DIR}/shared")

# surfker
add_subdirectory("${SOURCE_DIR}/surfker")

file(GLOB SOURCES "${SOURCE_DIR}/*.f90")

include_directories(
  ${PROJECT_SOURCE_DIR}/include
)

set(FILTERED_SOURCE_FILES "")
foreach(FILE ${SOURCES})
    if (NOT ${FILE} MATCHES "${SOURCE_DIR}/rsopt*.f90")
        list(APPEND FILTERED_SOURCE_FILES ${FILE})
    endif()
endforeach()

add_library(surfatt_modules OBJECT ${FILTERED_SOURCE_FILES})
target_link_libraries(surfatt_modules PRIVATE fortran-yaml shared_obj surfker_obj)
