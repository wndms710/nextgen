include ../../Makefile.in
-include $(CANTERA_HOME)/include/cantera/Cantera.mak

## check if cantera version < 2.3

ifdef CANTERA_VERSION
CT_VER_MAJOR := $(shell echo $(CANTERA_VERSION)| cut -f1 -d.)
CT_VER_MINOR := $(shell echo $(CANTERA_VERSION)| cut -f2 -d.)
CT_LT_2_3 := $(shell [ $(CT_VER_MAJOR) -lt 2 -o \( $(CT_VER_MAJOR) -eq 2 -a $(CT_VER_MINOR) -lt 3 \) ] && echo true)
ifeq ($(CT_LT_2_3),true)
  CXXFLAGS_NO_MPI += -DOLD_CANTERA
endif
CXXFLAGS_NO_MPI += -std=c++11
endif

NVCCFLAGS += -I../core/include
CXXFLAGS  += -I../core/include
CXXFLAGS  += -I../core/hdrlib
CXXFLAGS  += -I../core/common
CXXFLAGS  += -I../core/encrypt
CTI_LIBS = -L../core/common -lcti_core_s -L../core/encrypt -ltomlight  

CXXFLAGS_NO_MPI += $(CANTERA_INCLUDES)
CXXFLAGS_NO_MPI += -I../core/include
CXXFLAGS_NO_MPI += -I../core/hdrlib
CXXFLAGS_NO_MPI += -I../core/common
CXXFLAGS_NO_MPI += -I../core/encrypt

CXXFLAGS_NO_MPI += -O3 -g
CXXFLAGS += -O3 -g

default:		libchemistry.a createChemtable

libchemistry.a:		Chemtable.o
	@echo
	@echo "-------------------------------------------------"
	@echo "making chemistry components for charles..."
	@echo "-------------------------------------------------"
	@echo
	$(AR) $@ $^
	$(RANLIB) $@

createChemtable:
ifndef CANTERA_VERSION
	@echo "----------------------------------------------------------------"
	@echo " -- cantera not found , unable to build createChemtable" 
	@echo "----------------------------------------------------------------" && false
endif
	@echo
	@echo "-------------------------------------------------"
	@echo "making createChemtable..."
	@echo "CANTERA_VERSION = "$(CANTERA_VERSION)
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight.a
	make createChemtable.exe
	make ../../bin/createChemtable.exe


Chemtable.o:        Chemtable.hpp Chemtable.cpp
ChemtableGPU.o:     ChemtableGPU.hpp ChemtableGPU.cu
createChemtable.o:  ../core/lib/libcti_core_s.a ../core/encrypt/libtomlight.a Flamelet.hpp properties.hpp createChemtable.hpp createChemtable.cpp candle.hpp candle.o
PdfConvolution.o:   PdfConvolution.hpp PdfConvolution.cpp
candle.o:           candle.cpp candle.hpp StFlowLe.o Sim1DLe.o OneDimLe.o Domain1DLe.o refineLe.o boundaries1DLe.o MultiNewtonLe.o MultiJacLe.o
StFlowLe.o:         StFlowLe.cpp
Sim1DLe.o:          Sim1DLe.cpp
OneDimLe.o:         OneDimLe.cpp MultiNewtonLe.cpp
Domain1DLe.o:       Domain1DLe.cpp
refineLe.o:         refineLe.cpp
boundaries1DLe.o:   boundaries1DLe.cpp
MultiJacLe.o:       MultiJacLe.cpp
MultiNewtonLe.o:    MultiNewtonLe.cpp


libchemistry_gpu.a: Chemtable.o ChemtableGPU.o
	$(AR) $@ $^
	$(RANLIB) $@

candle: candle.cpp candle.o StFlowLe.o Sim1DLe.o OneDimLe.o Domain1DLe.o refineLe.o boundaries1DLe.o MultiNewtonLe.o MultiJacLe.o
	$(CXX_NO_MPI) -o candle.exe candle.o StFlowLe.o Sim1DLe.o OneDimLe.o Domain1DLe.o refineLe.o boundaries1DLe.o MultiNewtonLe.o MultiJacLe.o $(CANTERA_LIBS)

createChemtable.exe: ../core/common/libcti_core_s.a PdfConvolution_s.o createChemtable_s.o StFlowLe_s.o Sim1DLe_s.o OneDimLe_s.o Domain1DLe_s.o refineLe_s.o boundaries1DLe_s.o MultiNewtonLe_s.o MultiJacLe_s.o 
	$(CXX_NO_MPI) -o createChemtable.exe $^  $(CTI_LIBS) $(CLIBS_NO_MPI) $(CANTERA_LIBS)

../../bin/%.exe : %.exe
	cp $^ ../../bin/


%_s.o : %.cpp
	$(Q)$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI $< -c -o $@
	@echo "     CXX_S $<"

%_s.exe : %_s.o
	$(CXX_NO_MPI) -o $@ $^ -L../core/lib -lcti_core_s $(CLIBS_NO_MPI) $(CANTERA_LIBS)

#===============
# gpu lookup test
#===============
test.o: 		test.cu
	$(Q)$(NVCC) $(NVCCFLAGS) -c -o $@ $<
	
test.exe:		test.o Chemtable.o ChemtableGPU.o
	$(CXX) -o $@ $^ -L../core/lib -lcti_core -L$(CUDA_HOME)/lib64 -lcudart

clean:
	rm -vf *.o
	rm -vf *.a
	rm -vf *.exe
