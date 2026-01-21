# ===================================================
# Build a custom version of charles or postpro. Copy
# this file to your run directory and rename it 
# simply "Makefile".  Set an environment variable 
# named CTI_HOME pointing to the latest 
# release directory or specify the path below:
# 
CTI_HOME=../../../../
#
# ===================================================
include $(CTI_HOME)/Makefile.in

CXXFLAGS +=	-I$(CTI_HOME)/src/core/common
CXXFLAGS +=	-I$(CTI_HOME)/src/core/encrypt
CXXFLAGS +=	-I$(CTI_HOME)/src/surfer
CXXFLAGS +=	-I$(CTI_HOME)/src/stitch
CXXFLAGS += 	-I$(CTI_HOME)/src/chemistry
CXXFLAGS += 	-I$(CTI_HOME)/src/charles
CXXFLAGS += 	-I$(CTI_HOME)/src/charles/ideal_gas
CXXFLAGS += 	-I$(CTI_HOME)/src/charles/helmholtz
CXXFLAGS += 	-I$(CTI_HOME)/src/charles/premixed
CXXFLAGS += 	-I$(CTI_HOME)/src/charles/nonpremixed
CXXFLAGS += 	-I$(CTI_HOME)/src/charles/lsp
CXXFLAGS +=	-I$(CTI_HOME)/src/quiver
CXXFLAGS +=	-I$(CTI_HOME)/src/vof

CHARLES_OPT_FLAGS 	?=	
CXXFLAGS 	  	+=	$(CHARLES_OPT_FLAGS)

# by default, we'll use the charles opt flags for the out 
# of directory build
CXXFLAGS +=	$(CHARLES_OPT_FLAGS)

# this should also be updated when the core libraries are renamed..
CTI_LIBS = 	-L$(CTI_HOME)/src/core/common -lcti_core \
		-L$(CTI_HOME)/src/core/encrypt -ltomlight \
		-L$(CTI_HOME)/src/charles/lsp -lcti_lsp $(CLIBS)

CXXFLAGS_NO_MPI +=	-I$(CTI_HOME)/src/core/common
CXXFLAGS_NO_MPI +=	-I$(CTI_HOME)/src/core/encrypt
CXXFLAGS_NO_MPI +=	-I$(CTI_HOME)/src/post

CTI_LIBS_NO_MPI =	-L$(CTI_HOME)/src/core/common -lcti_core_s \
			-L$(CTI_HOME)/src/core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

charles:
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles.exe..."
	@echo "-------------------------------------------------"
	@echo
	make core_libs
	make charles.exe

core_libs:
	make -C $(CTI_HOME)/src/core/common libcti_core.a 
	make -C $(CTI_HOME)/src/core/encrypt

lsp_libs:
	make -C $(CTI_HOME)/src/charles/lsp cti_lsp

charles.cpp:
	cp $(CTI_HOME)/src/charles/charles.cpp .

$(CTI_HOME)/src/charles/%.o:
	#@echo $(CXX) $(subst $(CTI_HOME)/src/charles/,,$@)
	make -C $(CTI_HOME)/src/charles/ $(subst $(CTI_HOME)/src/charles/,,$@) 

CORE_HDR_DEPS=			$(CTI_HOME)/src/core/common/StaticSolver.hpp

SOLVER_HDR_DEPS=		$(CTI_HOME)/src/charles/States.hpp \
				$(CTI_HOME)/src/charles/FlowSolver.hpp \
				$(CTI_HOME)/src/charles/CtiMemory.hpp \
				$(CTI_HOME)/src/charles/CommContainers.hpp \
				$(CTI_HOME)/src/charles/Cht.hpp \
				$(CTI_HOME)/src/charles/BoundaryLayerDataExchanger.hpp

IG_HDR_DEPS=			$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolver.hpp \
				$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolverFlux.hpp \
				$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolverBcs.hpp \
				$(CTI_HOME)/src/charles/it/InflowTurbulenceIG.hpp 

PREMIXED_HDR_DEPS=		$(CTI_HOME)/src/chemistry/Chemtable.hpp \
				$(CTI_HOME)/src/charles/premixed/PremixedSolver.hpp \
				$(CTI_HOME)/src/charles/premixed/PremixedSolverFlux.hpp \
				$(CTI_HOME)/src/charles/premixed/PremixedSolverBcs.hpp \
				$(CTI_HOME)/src/charles/it/InflowTurbulenceP.hpp 

NONPREMIXED_HDR_DEPS=		$(CTI_HOME)/src/chemistry/Chemtable.hpp \
				$(CTI_HOME)/src/charles/nonpremixed/NonpremixedSolver.hpp \
				$(CTI_HOME)/src/charles/premixed/PremixedSolverFlux.hpp \
				$(CTI_HOME)/src/charles/nonpremixed/NonpremixedSolverBcs.hpp

HELMHOLTZ_HDR_DEPS=		$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolver.hpp \
				$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolverBcs.hpp 

IG_OBJECTS=			$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolver.o \
				$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolverSgs.o \
				$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolverIo.o \
				$(CTI_HOME)/src/charles/ideal_gas/IdealGasSolverBcs.o

PREMIXED_OBJECTS=		$(CTI_HOME)/src/charles/premixed/PremixedSolver.o \
				$(CTI_HOME)/src/charles/premixed/PremixedSolverSgs.o \
				$(CTI_HOME)/src/charles/premixed/PremixedSolverIo.o \
				$(CTI_HOME)/src/charles/premixed/PremixedSolverBcs.o

NONPREMIXED_OBJECTS=		$(CTI_HOME)/src/charles/nonpremixed/NonpremixedSolver.o \
				$(CTI_HOME)/src/charles/nonpremixed/NonpremixedSolverSgs.o \
				$(CTI_HOME)/src/charles/nonpremixed/NonpremixedSolverIo.o \
				$(CTI_HOME)/src/charles/nonpremixed/NonpremixedSolverBcs.o

HELMHOLTZ_OBJECTS=		$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolver.o \
				$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolverSgs.o \
				$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolverWm.o \
				$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolverIo.o \
				$(CTI_HOME)/src/charles/helmholtz/HelmholtzSolverBcs.o

charles.o :			$(CORE_HDR_DEPS) $(SOLVER_HDR_DEPS) $(IG_HDR_DEPS) \
				$(HELMHOLTZ_HDR_DEPS) $(PREMIXED_HDR_DEPS) $(NONPREMIXED_HDR_DEPS) \
				charles.cpp

CHARLES_OBJECTS = 		charles.o $(IG_OBJECTS) $(HELMHOLTZ_OBJECTS) \
				$(PREMIXED_OBJECTS) $(NONPREMIXED_OBJECTS)

charles.exe : 			$(CHARLES_OBJECTS) 
	make core_libs
	make lsp_libs
	$(CXX) -o $@ $^ $(CTI_LIBS)

# ----------------------------------------------------

clean:
	rm -vf charles.o 
	rm -vf charles.exe 

.PHONY: charles core_libs lsp_libs
