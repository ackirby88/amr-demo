CMAKE_MINIMUM_REQUIRED(VERSION 3.6)
PROJECT(amr-demo)

enable_language(Fortran)

# --------------------------------------------------------------------------- #
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
# --------------------------------------------------------------------------- #
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# --------------------------------------------------------------------------- #
# set C++17 standard for OCCA compatibility
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# --------------------------------------------------------------------------- #
# force GNU linker to include all shared libraries
if (NOT APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed")
else ()
  set(CMAKE_MACOSX_RPATH ON)
endif ()

# --------------------------------------------------------------------------- #
set(p4est_dir "" CACHE PATH "p4est install dir")
set(metis_dir "" CACHE PATH "metis install dir")
# --------------------------------------------------------------------------- #

#add_definitions(-Wall)
# --------------------------------------------------------------------------- #

add_subdirectory(solver/src)
add_subdirectory(physics/src)
add_subdirectory(amr/src)
