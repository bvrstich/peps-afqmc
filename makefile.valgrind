###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = peps
CPPSRC	= main.cpp\
           Random.cpp\
           Lattice.cpp\
           PEPS.cpp\
           MPS.cpp\
           MPO.cpp\
           compress.cpp\
           Global.cpp\
           Environment.cpp\
           Heisenberg.cpp\
           Trotter.cpp\
           propagate.cpp\
           btas_defs.cpp

OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

BTASLIB= /home/bright/btas/lib

INCLUDE = ./include

#LIBS= -lpthread -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBS= -lblas -llapacke

CC	= gcc
CXX	= g++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= -I$(INCLUDE) -g -std=c++11 -D_HAS_CBLAS -D_HAS_LAPACKE
LDFLAGS	= -g

# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME)
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@

# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

# -----------------------------------------------------------------------------
#   Make new documentation using doxygen
# -----------------------------------------------------------------------------
doc:
	@doxygen doc-config

# ====================== End of file 'makefile.in' ========================== #
