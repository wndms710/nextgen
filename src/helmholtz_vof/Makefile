include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir) 
CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt
CXXFLAGS +=	-I../charles
CXXFLAGS +=	-I../charles/lsp

CHARLES_OPT_FLAGS ?=	
CXXFLAGS 	  +=	$(CHARLES_OPT_FLAGS)

CTI_LIBS = 	-L../core/common -lcti_core -L../core/encrypt -ltomlight -L../charles/lsp -lcti_lsp $(CLIBS)

default:	charles_vof

../../bin/%.exe : %.exe
	cp $^ ../../bin/

charles_vof:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles_vof..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt 
	make -C ../charles/lsp 
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles_vof.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles_vof.exe
	make ../../bin/charles_vof.exe


CORE_HDR_DEPS=		../core/common/StaticSolver.hpp ../core/common/CtiScene.hpp
VOF_HDR_DEPS=		HelmholtzVofSolver.hpp HelmholtzVofSolverBcs.hpp LiquidFlux.hpp ../charles/FlowSolver.hpp

#---------------------------------
# boundary conditions
#---------------------------------
TrilinosInterface.o:            TrilinosInterface.hpp TrilinosInterface.cpp
HelmholtzVofSolverBcs.o:		$(CORE_HDR_DEPS) $(VOF_HDR_DEPS) HelmholtzVofSolverBcs.cpp 
LiquidFlux.o:           LiquidFlux.hpp LiquidFlux.cpp

#---------------------------------
# main solver driver
#---------------------------------

VOF_OBJECTS= 		charles_vof.o HelmholtzVofSolverBcs.o TrilinosInterface.o LiquidFlux.o

charles_vof.o :			../core/common/libcti_core.a $(CORE_HDR_DEPS) $(VOF_HDR_DEPS) charles_vof.cpp

charles_vof.exe :	$(VOF_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

clean:
	rm -vf *.o 
	rm -vf charles_vof.exe ../../bin/charles_vof.exe

.PHONY: charles_vof
