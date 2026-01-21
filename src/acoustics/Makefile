include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir) 
CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt

ACOUSTICS_OPT_FLAGS ?=	
CXXFLAGS +=	$(ACOUSTICS_OPT_FLAGS)

CTI_LIBS = 	-L../core/common -lcti_core -L../core/encrypt -ltomlight $(CLIBS)

default:	fwh resonance

fwh:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for fwh..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt 
	@echo
	@echo "-------------------------------------------------"
	@echo "making fwh.exe..."
	@echo "-------------------------------------------------"
	@echo
	make fwh.exe
	make ../../bin/fwh.exe

resonance: 
	@echo 
	@echo "-------------------------------------------------"
	@echo "making components for resonance..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt 
	@echo
	@echo "-------------------------------------------------"
	@echo "making resonance.exe..."
	@echo "-------------------------------------------------"
	@echo
	make resonance.exe


../../bin/%.exe : %.exe
	cp $^ ../../bin/

fwh.o :			fwh.cpp fwh.hpp observer.hpp

resonance.o :		resonance.cpp ../core/common/StaticSolver.hpp

%.exe :			%.o	
	$(CXX) -o $@ $^ $(CTI_LIBS)

clean:
	rm -vf *.o 
	rm -vf fwh.exe ../../bin/fwh.exe
	rm -vf resonance.exe ../../bin/resonance.exe

.PHONY: fwh resonance
