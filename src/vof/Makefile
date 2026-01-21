include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir) 
CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt
CXXFLAGS +=	-I../charles
CXXFLAGS +=	-I../charles/lsp

CHARLES_OPT_FLAGS ?=	
CXXFLAGS 	  +=	$(CHARLES_OPT_FLAGS)

CTI_LIBS = 	-L../core/common -lcti_core -L../core/encrypt -ltomlight -L../charles/lsp -lcti_lsp $(CLIBS)

default:	vof

../../bin/%.exe : %.exe
	cp $^ ../../bin/

vof:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for vof..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt 
	make -C ../charles/lsp 
	@echo
	@echo "-------------------------------------------------"
	@echo "making vof.exe..."
	@echo "-------------------------------------------------"
	@echo
	make vof.exe
	make ../../bin/vof.exe


CORE_HDR_DEPS=		../core/common/StaticSolver.hpp ../core/common/CtiScene.hpp
VOF_HDR_DEPS=		VofSolver.hpp VofSolverBcs.hpp DropStats.hpp ../charles/FlowSolver.hpp

#---------------------------------
# boundary conditions
#---------------------------------
TrilinosInterface.o:            TrilinosInterface.hpp TrilinosInterface.cpp
VofSolverBcs.o:		$(CORE_HDR_DEPS) $(VOF_HDR_DEPS) VofSolverBcs.cpp 
DropStats.o:            DropStats.hpp DropStats.cpp

#---------------------------------
# main solver driver
#---------------------------------

VOF_OBJECTS= 		vof.o VofSolverBcs.o TrilinosInterface.o DropStats.o

vof.o :			../core/common/libcti_core.a $(CORE_HDR_DEPS) $(VOF_HDR_DEPS) vof.cpp

vof.exe :	$(VOF_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

clean:
	rm -vf *.o 
	rm -vf vof.exe ../../bin/vof.exe

.PHONY: vof
