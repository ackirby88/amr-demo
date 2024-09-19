file(GLOB PHYS_H  *.h)
set(PHYS_SRC
    physics_mpi.c
    physics_runtime_statistics.c
    physics_utilities.c
    physics_norm_residual.c
    physics_initialize.c
    physics_refine_operator_2d.c
    physics_refine_operator_3d.c
    physics_coarsen_operator_2d.c
    physics_coarsen_operator_3d.c
    physics_residual_subgrid.c
    physics_residual_bc_unstructured.c
    physics_residual_bc_structured.c
    physics_residual_face_full_unstructured.c
    physics_residual_face_full_structured.c
    physics_residual_face_hang_unstructured.c
    physics_residual_face_full_hang_structured.c
    physics_residual_face_hang_full_structured.c
    physics_rhs_unstructured.c
    physics_rhs_structured.c
    physics_rhs.c
    physics_time_step.c
    rk_time_integrators.c
    physics_runge_kutta_interface.c
    physics_solve_interface.c
    physics_go.c
    physics_tag_feature.c
    physics_regrid_points.c
    physics_evolve_interface.c
    physics_initialize_quadrant.c
    hpc_amr_initialize_solver.c
    hpc_amr_initialize_quadrant.c
    hpc_amr_refine_operator_2d.c
    hpc_amr_refine_operator_3d.c
    hpc_amr_coarsen_operator_2d.c
    hpc_amr_coarsen_operator_3d.c
    hpc_amr_tag_feature.c
    hpc_amr_regrid_points.c
    hpc_amr_evolve_solver.c
)

add_library(hpc_phys SHARED ${PHYS_SRC} ${PHYS_H})
target_link_libraries(hpc_phys hpc_solver)

install(TARGETS hpc_phys
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib/lantern 
        ARCHIVE DESTINATION lib/lantern
)

INSTALL(FILES ${PHYS_H} DESTINATION include/lantern/hpc)
