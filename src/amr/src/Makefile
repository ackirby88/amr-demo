CC      = mpicc
CFLAGS  = -O2 -c
#CFLAGS  = -c -g -O0 -Wall

#CFLAGS += -D_2D_

P4EST_DIR = /home/akirby/Software/p4est-2.0/build/local

LDFLAGS  = -L$(P4EST_DIR)/lib/ -Wl,-rpath=$(P4EST_DIR)/lib/ -lp4est
LDFLAGS += -L$(P4EST_DIR)/lib/ -Wl,-rpath=$(P4EST_DIR)/lib/ -lsc
LDFLAGS += -lm -lrt

INCLUDES= -I$(P4EST_DIR)/include


LIBMETIS = #-L/usr/local/lib -Wl,-rpath=/usr/local/lib -lmetis


PHYSICS_PATH = /home/akirby/bin
PHYSICS_LIB  = physics
SOLVER_LIB   = fvsolver

LIBSOLVER   = -L$(PHYSICS_PATH) -Wl,-rpath=$(PHYSICS_PATH) -l$(PHYSICS_LIB) -l$(SOLVER_LIB) -lstdc++

SOURCES =   \
    amr_utilities.c \
    amr_p4est_utilities.c \
    amr_tag.c \
    amr_mpi.c \
    amr_vtk.c \
    amr_external.c \
    amr_initialize.c \
    amr_regrid.c \
    amr_evolve.c \
    main.c \

EXECUTABLE = hpc_driver


OBJECTS= $(addsuffix .o,$(basename $(SOURCES)))


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) $(LIBMETIS) $(LIBSOLVER) -o ~/bin/$(EXECUTABLE)


.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean:
	$(RM) /home/akirby/bin/$(EXECUTABLE).mpi *.o
