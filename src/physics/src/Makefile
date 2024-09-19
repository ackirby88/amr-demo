CC      = mpicc
CXX	= mpicxx
CFLAGS  = -O2 -c -fPIC
#CFLAGS  = -c -g -O0 -Wall -fPIC


LD	= $(CC)

LIBS = -L/home/akirby/bin -Wl,-rpath=/home/akirby/bin/ -lfvsolver


SOURCES =   \
    physics_mpi.c \
    physics_runtime_statistics.c \
    physics_utilities.c \
    physics_norm_residual.c \
    physics_initialize.c \
    physics_refine_operator_2d.c \
    physics_refine_operator_3d.c \
    physics_coarsen_operator_2d.c \
    physics_coarsen_operator_3d.c \
    physics_residual_subgrid.c \
    physics_residual_bc_unstructured.c \
    physics_residual_bc_structured.c \
    physics_residual_face_full_unstructured.c \
    physics_residual_face_full_structured.c \
    physics_residual_face_hang_unstructured.c \
    physics_residual_face_full_hang_structured.c \
    physics_residual_face_hang_full_structured.c \
    physics_rhs_unstructured.c \
    physics_rhs_structured.c \
    physics_rhs.c \
    physics_time_step.c \
    rk_time_integrators.c \
    physics_runge_kutta_interface.c \
    physics_solve_interface.c \
    physics_go.c \
    physics_tag_feature.c \
    physics_regrid_points.c \
    physics_evolve_interface.c \
    physics_initialize_quadrant.c \
    hpc_amr_initialize_solver.c \
    hpc_amr_initialize_quadrant.c \
    hpc_amr_refine_operator_2d.c \
    hpc_amr_refine_operator_3d.c \
    hpc_amr_coarsen_operator_2d.c \
    hpc_amr_coarsen_operator_3d.c \
    hpc_amr_tag_feature.c \
    hpc_amr_regrid_points.c \
    hpc_amr_evolve_solver.c \

EXECUTABLE = physics_main


OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))
DEFS = $(foreach define, $(DEFINE),-D$(define))


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(LD) $(OBJECTS) $(LIBS) $(LDCFLAGS) -shared -o ~/bin/libphysics.so -lc
		@echo ""
		@echo " >> Status : ok! "
		@echo ""


.SUFFIXES: .c .C .o

.c.o:
	$(CC) $(CFLAGS) $(DEFS) $< -o $@
	
.C.o:
	$(CXX) $(CFLAGS) $(DEFS) $< -o $@

clean:
	$(RM) ~/bin/libphysics.so *.o

