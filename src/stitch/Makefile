include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir)
CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt

STITCH_OPT_FLAGS ?=
CXXFLAGS +=	$(STITCH_OPT_FLAGS)

CTI_LIBS = 	-L../core/common -lcti_core -L../core/encrypt -ltomlight $(CLIBS)

default:	stitch

../../bin/%.exe : %.exe
	cp $^ ../../bin/

# ======================== stitch ===============================

stitch:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for stitch..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making stitch.exe..."
	@echo "-------------------------------------------------"
	@echo
	make stitch.exe
	make ../../bin/stitch.exe

DEP = ../core/common/libcti_core.a ../core/encrypt/libtomlight.a

PartData.o:	$(DEP) ../core/common/StZone.hpp ../core/common/SurfaceShm.hpp PartData.hpp Points.hpp PartData.cpp

StitchNS.o:	$(DEP) PartData.hpp HcpPointBuilder.hpp StitchNS.hpp StitchNS.cpp

CuttableVoronoiData.o:	$(DEP) CuttableVoronoiData.hpp CuttableVoronoiData.cpp
CuttableVoronoiData.o:	$(DEP) corner_flag1.hpp corner_flag2.hpp corner_flag3.hpp cvd_doit.hpp doit.hpp doit2.hpp ice_backward.hpp ice_forward.hpp \
	CuttableVoronoiData.hpp CuttableVoronoiData.cpp

stitch.o :	$(DEP) PartData.hpp StitchNS.hpp \
	compressAndSend.hpp HcpPointBuilder.hpp \
	stitch.siData stitch.cpp

stitch.exe :	PartData.o CuttableVoronoiData.o StitchNS.o stitch.o
	$(CXX) -o $@ $^ $(CTI_LIBS)

#========================= pbin generation utility =====================

generate_pbin.exe:	generate_pbin.cpp ../core/common/Adt.hpp
	$(CXX_NO_MPI) -o $@ $(CXXFLAGS_NO_MPI) -I../core/common $<

generate_pbin2.o: 	generate_pbin2.cpp ../core/common/Adt.hpp \
			../core/common/MpiStuff.hpp ../core/common/DistributedDataExchanger.hpp \
			../core/common/CommonIo.hpp ../core/common/Macros.hpp ../core/common/MiscUtils.hpp

generate_pbin2.exe:	generate_pbin2.o ../core/common/MiscUtils.o ../core/common/MpiStuff.o
	$(CXX) -o $@ $^ $(CLIBS)

generate_pbin2:
	make -C ../core/common MpiStuff.o
	make -C ../core/common MiscUtils.o
	make generate_pbin2.exe

generate_pbin3.exe:	generate_pbin3.o
	$(CXX) -o $@ $^ $(CLIBS)

generate_pbin3:
	make generate_pbin3.exe


#tetMshToTetBin.exe:	tetMshToTetBin.cpp ../core/common/FluentReader.hpp \
#	../core/common/NoMpi_s.o ../core/common/MpiStuff_s.o ../core/common/MiscUtils_s.o \
#	../core/common/Params_s.o ../core/common/CTI_s.o ../core/common/CtiRegister_s.o

#tetMshToTetBin:
#	make -C ../core/common NoMpi_s.o
#	make -C ../core/common MpiStuff_s.o
#	make -C ../core/common MiscUtils_s.o
#	make tetMshToTetBin.exe

tetMshToTetBin:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for tetMshToTetBin..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making tetMshToTetBin.exe..."
	@echo "-------------------------------------------------"
	@echo
	make tetMshToTetBin.exe
	make ../../bin/tetMshToTetBin.exe

tetMshToTetBin_s.o:	../core/common/FluentReader.hpp tetMshToTetBin.cpp

tetMshToTetBin.exe:	../core/common/libcti_core_s.a tetMshToTetBin.o
	$(CXX_NO_MPI) -o $@ $^ -L../core/common -lcti_core_s -L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)


fluentMshToPbin:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for fluentMshToPbin..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt libtomlight.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making fluentMshToPbin.exe..."
	@echo "-------------------------------------------------"
	@echo
	make fluentMshToPbin.exe
	make ../../bin/fluentMshToPbin.exe

fluentMshToPbin.o:	../core/common/DataContainer.hpp fluentMshToPbin.cpp

fluentMshToPbin.exe:	../core/common/libcti_core.a fluentMshToPbin.o
	$(CXX) -o $@ $^ -L../core/common -lcti_core -L../core/encrypt -ltomlight $(CLIBS)

#	$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI -o $@ $^ $(CLIBS_NO_MPI)




# ================================= tests ===============================

VoxelGrid.o: ../core/common/SurfaceShm.hpp VoxelGrid.hpp
testVoxelGrid.o: VoxelGrid.o testVoxelGrid.cpp
testVoxelGrid.exe: testVoxelGrid.o VoxelGrid.o
	$(CXX) -o $@ $^ $(CTI_LIBS)
testVoxelGrid:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for testVoxelGrid..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making testVoxelGrid.exe..."
	@echo "-------------------------------------------------"
	@echo
	make testVoxelGrid.exe

testHcp2d.o: ../core/common/SurfaceShm.hpp testHcp2d.cpp
testHcp2d.exe: testHcp2d.o
	$(CXX) -o $@ $^ $(CTI_LIBS)
testHcp2d:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for testHcp2d..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making testHcp2d.exe..."
	@echo "-------------------------------------------------"
	@echo
	make testHcp2d.exe


# ==================================------===============================
# serial stitch...

CXXFLAGS_NO_MPI +=	-I../core/common
CXXFLAGS_NO_MPI +=	-I../core/encrypt

CTI_LIBS_S = 	-L../core/common -lcti_core_s -L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

#Rule for serial compilation
%_s.o : %.cpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@
	@#@echo "     CXX_S $<"

#../core/common/libcti_core_s.a:
#	make -C ../core/common libcti_core_s.a

stitch_s:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for stitch_s..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making stitch.exe..."
	@echo "-------------------------------------------------"
	@echo
	make stitch_s.exe

DEP_S = ../core/common/libcti_core_s.a ../core/encrypt/libtomlight_s.a

PartData_s.o:	$(DEP_S) ../core/common/StZone.hpp ../core/common/SurfaceShm.hpp PartData.hpp Points.hpp PartData.cpp

StitchNS_s.o:	$(DEP_S) PartData.hpp StitchNS.hpp StitchNS.cpp

CuttableVoronoiData_s.o:	$(DEP_S) corner_flag1.hpp corner_flag2.hpp corner_flag3.hpp cvd_doit.hpp doit.hpp doit2.hpp ice_backward.hpp ice_forward.hpp \
	CuttableVoronoiData.hpp CuttableVoronoiData.cpp

stitch_s.o :	$(DEP_S) PartData.hpp StitchNS.hpp \
	compressAndSend.hpp HcpPointBuilder.hpp \
	stitch.siData stitch.cpp

stitch_s.exe :	PartData_s.o CuttableVoronoiData_s.o StitchNS_s.o stitch_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS_S)

# ==================================------===============================

cvd_debug:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for cvd_debug..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making cvd_debug.exe..."
	@echo "-------------------------------------------------"
	@echo
	make cvd_debug.exe

cvd_debug_s.o :	$(DEP_S) cvd_debug.cpp

cvd_debug.exe :	CuttableVoronoiData_s.o cvd_debug_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS_S)

# ==================================------===============================

tests:
	make -C tests

# ==================================------===============================

clean:
	rm -vf *.o
	rm -vf stitch.exe stitch_s.exe testVoxelGrid.exe testHcp2d.exe \
		generate_pbin.exe generate_pbin2.exe generate_pbin3.exe tetMshToTetBin.exe \
		cvd_debug.exe
	make -C tests clean

.PHONY: stitch generate_pbin2 generate_pbin3 testVoxelGrid testHcp2d tests tetMshToTetBin
