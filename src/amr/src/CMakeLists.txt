# ======================================= #
# Library header files to be included
include_directories(${p4est_dir}/include
                    ${metis_dir}/include
                    ./include)

# Libraries to be linked
find_library(p4est_lb     NAMES p4est           PATHS ${p4est_dir}/lib)
find_library(sc_lb        NAMES sc              PATHS ${p4est_dir}/lib)
find_library(metis_lb     NAMES metis           PATHS ${metis_dir}/lib)
# ======================================= #

# Uncomment for 2D simulations
add_definitions(-D_2D_)

file(GLOB AMR_H  *.h)
set(AMR_SRC 
    amr_utilities.c
    amr_p4est_utilities.c
    amr_tag.c
    amr_mpi.c
    amr_vtk.c
    amr_external.c
    amr_initialize.c
    amr_regrid.c
    amr_evolve.c
)

add_library(hpc_amr SHARED ${AMR_SRC} ${AMR_H})
target_link_libraries(hpc_amr ${sc_lb} ${p4est_lb} ${metis_lb} ${occa_lb} m)


add_executable(demo.mpi main.c)
target_link_libraries(demo.mpi hpc_amr hpc_phys hpc_solver ${sc_lb} ${p4est_lb} ${metis_lb} ${occa_lb} m)

install(TARGETS hpc_amr demo.mpi
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib/ 
        ARCHIVE DESTINATION lib/
)
