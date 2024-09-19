FC      = mpif90
FFLAGS  = -O2 -c -fPIC
#FFLAGS  = -c -g -O0 -Wall -fPIC -fbacktrace

LD = $(FC)

SOURCES =   \
    my_kinddefs.F90 \
    cell_volume_module.F90 \
    tagging_methods_module.F90 \
    subgrid_index_2d_module.F90 \
    subgrid_index_3d_module.F90 \
    bilinear_node_2d_module.F90 \
    bilinear_node_3d_module.F90 \
    integrate_2d_module.F90 \
    integrate_3d_module.F90 \
    normals_2d_module.F90 \
    normals_3d_module.F90 \
    flux_2d_module.F90 \
    flux_3d_module.F90 \
    bc_2d_module.F90 \
    bc_3d_module.F90 \
    initial_boundary_2d_module.F90 \
    initial_boundary_3d_module.F90 \
    initial_condition_2d_module.F90 \
    initial_condition_3d_module.F90 \
    solver_time_step_2d_module.F90 \
    solver_time_step_3d_module.F90 \
    solver_time_step.F90 \
    solver_tag_feature.F90 \
    solver_initialize_quadrant_2d_module.F90 \
    solver_initialize_quadrant_3d_module.F90 \
    solver_initialize_quadrant.F90 \
    solver_refine_operator_2d.F90 \
    solver_refine_operator_3d.F90 \
    solver_coarsen_operator_2d.F90 \
    solver_coarsen_operator_3d.F90 \
    solver_residual_bc_unstructured_2d_module.F90 \
    solver_residual_bc_unstructured_3d_module.F90 \
    solver_residual_bc_unstructured.F90 \
    solver_residual_bc_structured_2d_module.F90 \
    solver_residual_bc_structured_3d_module.F90 \
    solver_residual_bc_structured.F90 \
    solver_residual_face_full_unstructured_2d_module.F90 \
    solver_residual_face_full_unstructured_3d_module.F90 \
    solver_residual_face_full_unstructured.F90 \
    solver_residual_face_full_structured_2d_module.F90 \
    solver_residual_face_full_structured_3d_module.F90 \
    solver_residual_face_full_structured.F90 \
    solver_residual_face_hang_unstructured_2d.F90 \
    solver_residual_face_hang_unstructured_3d.F90 \
    solver_residual_face_full_hang_structured_2d.F90 \
    solver_residual_face_full_hang_structured_3d.F90 \
    solver_residual_face_hang_full_structured_2d.F90 \
    solver_residual_face_hang_full_structured_3d.F90 \
    solver_residual_subgrid_2d_module.F90 \
    solver_residual_subgrid_3d_module.F90 \
    solver_residual_subgrid.F90 \

EXECUTABLE = amr_main 


OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))
DEFS = $(foreach define, $(DEFINE),-D$(define))


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(LD) $(OBJECTS) $(LIBS) $(LDCFLAGS) -shared -o ~/bin/libfvsolver.so -lc
		@echo ""
		@echo " >> Status : ok! "
		@echo ""


.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(FFLAGS) $(DEFS) $< -o $@
	
clean:
	$(RM) ~/bin/*.so *.o *.mod
