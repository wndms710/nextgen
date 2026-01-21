include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir) 
CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt

CTI_LIBS = 	-L../core/common -lcti_core \
		-L../core/encrypt -ltomlight $(CLIBS)

default:	potential

potential:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for potential..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt 
	@echo
	@echo "-------------------------------------------------"
	@echo "making potential.exe..."
	@echo "-------------------------------------------------"
	@echo
	make potential.exe

CORE_HDR_DEPS=	../core/common/StaticSolver.hpp ../core/common/StripedMesh.hpp

potential.o :	../core/common/libcti_core.a $(CORE_HDR_DEPS) PotentialSolver.hpp potential.cpp

potential.exe :	potential.o
	$(CXX) -o $@ $^ $(CTI_LIBS)

clean:
	rm -vf *.o 
	rm -vf potential.exe

.PHONY: potential
