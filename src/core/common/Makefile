include ../../../Makefile.in
include ../../../makefiles/Makefile.version

CXXFLAGS 	+= -I../encrypt
CXXFLAGS 	+= -I.

CORE_OPT_FLAGS ?=
CXXFLAGS        += $(CORE_OPT_FLAGS)

BLAS_LIB ?= -lblas
LAPACK_LIB ?= -llapack

CTI_VERSION ?= git describe --all --long --always
DOCS_VERSION ?= master

#
CTI_DEFS=	Defs.hpp ByteSwap.hpp

default:	libcti_core.a


####################################################################################################
# parallel and serial core libraries...

CORE_OBJS=	NoMpi.o MpiStuff.o Guardian.o MiscUtils.o Params.o CtiRegister.o SimpleFunc.o CTI.o \
	Version.o AverageOperator.o RestartHashUtilities.o Logger.o PeriodicData.o CtiCanvas.o LesImageMapper.o ColorMap.o SimpleGeom.o \
	Stl.o Multigrid.o StaticSolver.o StaticSolver_multigrid.o StaticSolver_probes.o StaticSolver_write_data.o WebUI.o \
	GeomUtils.o SplineStuff.o GeodesicSphere.o 

libcti_core.a: $(CORE_OBJS)
	$(AR) $@ $^
	$(RANLIB) $@

CORE_OBJS_S=	NoMpi_s.o MpiStuff_s.o Guardian_s.o MiscUtils_s.o Params_s.o CtiRegister_s.o CTI_s.o \
	Version_s.o RestartHashUtilities_s.o Logger_s.o PeriodicData_s.o CtiCanvas_s.o LesImageMapper_s.o ColorMap_s.o WebUI_s.o \
	GeomUtils_s.o SplineStuff_s.o GeodesicSphere_s.o SimpleGeom_s.o SimpleFunc_s.o 

libcti_core_s.a: $(CORE_OBJS_S)
	$(AR) $@ $^
	$(RANLIB) $@

####################################################################################################
# basic utilities
NoMpi.o:		$(CTI_DEFS) NoMpi.hpp NoMpi.cpp
MpiStuff.o: 		$(CTI_DEFS) NoMpi.hpp MpiStuff.hpp MpiStuff.cpp
Guardian.o: 		$(CTI_DEFS) Guardian.hpp Guardian.cpp
MiscUtils.o:		$(CTI_DEFS) MpiStuff.hpp MiscUtils.hpp MiscUtils.cpp
CtiRegister.o:		$(CTI_DEFS) MpiStuff.hpp CtiRegister.hpp CtiRegister.cpp
Params.o: 		$(CTI_DEFS) MpiStuff.hpp CtiRegister.hpp Params.hpp Params.cpp
SimpleFunc.o:		$(CTI_DEFS) Params.hpp SimpleFunc.hpp SimpleFunc.cpp
CTI.o:			$(CTI_DEFS) MpiStuff.o Guardian.o CtiRegister.o CTI.hpp CTI.cpp
AverageOperator.o:	$(CTI_DEFS) DistributedDataExchanger.hpp MpiStuff.o AverageOperator.hpp AverageOperator.cpp
Units.o:		$(CTI_DEFS) MpiStuff.hpp Units.hpp Units.cpp
RestartHashUtilities.o: $(CTI_DEFS) MpiStuff.hpp RestartHashUtilities.hpp RestartHashUtilities.cpp
Logger.o:               $(CTI_DEFS) Logger.hpp Logger.cpp
PeriodicData.o:		$(CTI_DEFS) PeriodicTransform.hpp PeriodicData.hpp PeriodicData.cpp
SimpleGeom.o:		CTI.hpp SimpleGeom.hpp SimpleGeom.cpp
Stl.o:			Stl.hpp Stl.cpp
WebUI.o:		$(CTI_DEFS) WebUIOutput.hpp WebUI.hpp WebUI.cpp
Multigrid.o:		$(CTI_DEFS) DistributedDataExchanger.hpp Utils.hpp MpiStuff.o MiscUtils.o
# geometry...
GeomUtils.o:		$(CTI_DEFS) Common.hpp SimpleTri.hpp GeomUtils.hpp GeomUtils.cpp
SplineStuff.o:		$(CTI_DEFS) Common.hpp SplineStuff.hpp SplineStuff.cpp
GeodesicSphere.o:	$(CTI_DEFS) Common.hpp GeodesicSphere.hpp GeodesicSphere.cpp

# dynamically grab the core version of the git repository so that
# it can be outputed ...
Version.o: ../../.git/HEAD CTI.hpp ../../../makefiles/Makefile.version
	@echo '#include "CTI.hpp"' > Version.cpp
	@echo 'const char* CTI::cti_core_version = "' | tr -d "\r\n" >> Version.cpp
	@$(CTI_VERSION) | tr -d "\r\n" >> Version.cpp
	@echo '";' >> Version.cpp
	@echo 'const char* CTI::cti_core_date = "' | tr -d "\r\n" >> Version.cpp
	@date | tr -d "\r\n" >> Version.cpp
	@echo '";' >> Version.cpp
	@echo 'const char* CTI::cti_docs_version = "$(DOCS_VERSION)";' >> Version.cpp
	$(CXX) $(CXXFLAGS)   -c -o $@ Version.cpp
	@#@echo   "     CXX Version.cpp"

# if this is not a git repo, then this creates .git/HEAD as an empty file...
../../.git/HEAD:
	mkdir ../../.git
	touch ../../.git/HEAD

##########################################################
#Serial Files
CXXFLAGS_NO_MPI  += -I./
CXXFLAGS_NO_MPI  += -I../encrypt

%_s.o : %.cpp
	$(Q)$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) $< -c -o $@
	@echo "     CXX_S $<"

#Rule for serial compilation
%_s.o : %.cpp
	$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI $< -c -o $@

#Needed for LicenseUtility
Guardian_s.o:		$(CTI_DEFS) Guardian.hpp Guardian.cpp
#Needed for SimulationTimeline, InputFileAssistant
Params_s.o:		$(CTI_DEFS) Params.hpp Params.cpp
RestartHashUtilities_s.o:	$(CTI_DEFS) MpiStuff.hpp RestartHashUtilities.hpp RestartHashUtilities.cpp
MpiStuff_s.o:		$(CTI_DEFS) NoMpi.hpp MpiStuff.hpp MpiStuff.cpp
NoMpi_s.o:		$(CTI_DEFS) NoMpi.hpp NoMpi.cpp
CTI_s.o:		$(CTI_DEFS) MpiStuff_s.o Guardian_s.o CTI.hpp CTI.cpp
MiscUtils_s.o:		$(CTI_DEFS) MpiStuff.hpp MiscUtils.hpp MiscUtils.cpp
Logger_s.o:            	$(CTI_DEFS) Logger.hpp Logger.cpp
PeriodicData_s.o:	$(CTI_DEFS) PeriodicTransform.hpp PeriodicData.hpp PeriodicData.cpp
CtiRegister_s.o:	$(CTI_DEFS) MpiStuff.hpp CtiRegister.hpp CtiRegister.cpp
SimpleGeom_s.o:		CTI.hpp SimpleGeom.hpp SimpleGeom.cpp
Stl_s.o:	 	Stl.hpp Stl.cpp
WebUI_s.o:		$(CTI_DEFS) WebUIOutput.hpp WebUI.hpp WebUI.cpp
Multigrid_s.o:		$(CTI_DEFS) DistributedDataExchanger.hpp Utils.hpp MpiStuff_s.o MiscUtils_s.o
# geometry...
GeomUtils_s.o:		$(CTI_DEFS) Common.hpp SimpleTri.hpp GeomUtils.hpp GeomUtils.cpp
SplineStuff_s.o:	$(CTI_DEFS) Common.hpp SplineStuff.hpp SplineStuff.cpp
GeodesicSphere_s.o:	$(CTI_DEFS) Common.hpp GeodesicSphere.hpp GeodesicSphere.cpp

#
Version_s.o: ../../.git/HEAD CTI.hpp ../../../makefiles/Makefile.version
	@echo '#include "CTI.hpp"' > Version.cpp
	@echo 'const char* CTI::cti_core_version = "' | tr -d "\r\n" >> Version.cpp
	@$(CTI_VERSION) | tr -d "\r\n" >> Version.cpp
	@echo '";' >> Version.cpp
	@echo 'const char* CTI::cti_core_date = "' | tr -d "\r\n" >> Version.cpp
	@date | tr -d "\r\n" >> Version.cpp
	@echo '";' >> Version.cpp
	@echo 'const char* CTI::cti_docs_version = "$(DOCS_VERSION)";' >> Version.cpp
	$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI  -c -o $@ Version.cpp
	@#@echo   "     CXX_NO_MPI Version.cpp"


#########################################################
#Testing

# --------------------------------------------
# the Cti environment including param parsing is
# tested with this routine...
# --------------------------------------------

testCti:
	make testCti.exe

testCti.o:	testCti.cpp
testCti.exe:	NoMpi.o MpiStuff.o Guardian.o MiscUtils.o Params.o CtiRegister.o CTI.o Version.o RestartHashUtilities.o \
	Logger.o PeriodicData.o ../encrypt/libtomlight.a testCti.o
	$(CXX) -o $@ $^ $(CLIBS)

# --------------------------------------------
# the StripedMesh class reads the restart (v5+)
# into memory striped...
# --------------------------------------------

testStripedMesh:
	make testStripedMesh.exe

testStripedMesh.o:	StripedMesh.hpp testStripedMesh.cpp
testStripedMesh.exe:	NoMpi.o MpiStuff.o Guardian.o MiscUtils.o Params.o CtiRegister.o CTI.o Version.o RestartHashUtilities.o \
	Logger.o PeriodicData.o ../encrypt/libtomlight.a testStripedMesh.o
	$(CXX) -o $@ $^ $(CLIBS)

# --------------------------------------------
# the StaticSolver is the base class for all
# static (i.e. non-moving) solvers...
# --------------------------------------------

testStaticSolver:
	make testStaticSolver.exe

PROBE_HEADERS = 	PointProbe.hpp MultiFluxProbe.hpp FluxProbe.hpp SimpleGeom.hpp ConditionalProbe.hpp \
			PointCloudProbe.hpp VolumetricProbe.hpp FazoneProbe.hpp PdfProbe.hpp FwhSurface.hpp
COMM_HEADERS  =		Prcomm.hpp updateCvDN.hpp updateCvDN3.hpp updateCvDN33.hpp updateCvDataSeparateGhosts.hpp \
			updateCv2DN.hpp updateCv2DN3.hpp updateCv2DN33.hpp updateCv2DNReverse.hpp updateNoData.hpp \
			updateFaDN.hpp updateFaDN3.hpp updateCcIData.hpp updateCcIDataReverse.hpp \
			updateCcData.hpp updateCbIData.hpp updateCfIData.hpp DataExchanger.hpp \
			updateSubSurfaceDN.hpp updateSubSurfaceDN3.hpp updateCvdDN.hpp updateCvdDN3.hpp
DATA_HEADERS  = 	DataContainer.hpp CasDatReader.hpp DataWriter.hpp StripedMesh.hpp SnapshotIO.hpp
VIZ_HEADERS   = 	CtiScene.hpp SimpleTri.hpp

SS_HEADERS    = 	StaticSolver.hpp KillfileReader.hpp CtiRegister.hpp SubSurface.hpp SurfaceShm.hpp tetCutPlane.hpp \
			tetCutPlaneData.hpp tetCutPlaneWgts.hpp CtiInterpolateTransform.hpp $(PROBE_HEADERS) $(COMM_HEADERS) \
			$(DATA_HEADERS) $(VIZ_HEADERS) FluentReader.hpp Octree.hpp

testStaticSolver.o:	libcti_core.a $(SS_HEADERS) testStaticSolver.cpp
testStaticSolver.exe:	testStaticSolver.o libcti_core.a ../encrypt/libtomlight.a
	$(CXX) -o $@ $^ $(CLIBS)

# --------------------------------------------
# the SurfaceShm is a parallel surface class...
# --------------------------------------------

testSurfaceShm:
	make testSurfaceShm.exe

testSurfaceShm.o:	SurfaceShm.hpp testSurfaceShm.cpp
testSurfaceShm.exe:	NoMpi.o MpiStuff.o Guardian.o MiscUtils.o Params.o CtiRegister.o CTI.o Version.o RestartHashUtilities.o \
	Logger.o PeriodicData.o ../encrypt/libtomlight.a testSurfaceShm.o
	$(CXX) -o $@ $^ $(CLIBS)

# --------------------------------------------

testExchanger:
	make testExchanger.exe

testExchanger.o:	$(SS_HEADERS) testExchanger.cpp
testExchanger.exe:	testExchanger.o libcti_core.a ../encrypt/libtomlight.a
	$(CXX) -o $@ $^ $(CLIBS)

testMesh:
	make testMesh.exe

testMesh.o:	libcti_core.a StripedMesh.hpp updateCvDN.hpp updateCvDN3.hpp updateCv2DN.hpp updateCv2DN3.hpp updateNoData.hpp DataContainer.hpp CasDatReader.hpp StaticSolver.hpp DataWriter.hpp CtiScene.hpp testMesh.cpp
testMesh.exe:	testMesh.o libcti_core.a ../encrypt/libtomlight.a
	$(CXX) -o $@ $^ $(CLIBS)

testRegister:
	make testRegister.exe

testRegister.o:		$(SS_HEADERS) testRegister.cpp
testRegister.exe:	testRegister.o libcti_core.a ../encrypt/libtomlight.a
	$(CXX) -o $@ $^ $(CLIBS)

testIO:
	make testIO.exe

testIO.o:		$(SS_HEADERS) testIO.cpp
testIO.exe:	testIO.o libcti_core.a ../encrypt/libtomlight.a
	$(CXX) -o $@ $^ $(CLIBS)


# --------------------------------------------
# test
#---------------------------------------------

test:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for test_core..."
	@echo "-------------------------------------------------"
	make -C .  libcti_core.a
	make -C ../encrypt libtomlight.a
	@echo
	make test_core.exe

tests/testParams.o: tests/testParams.cpp
tests/testStripedMesh.o: tests/testStripedMesh.cpp

TEST_OBJS= tests/testStripedMesh.o \
           tests/testParams.o

tests/test_core.o: ../catch1/catch.hpp tests/test_core.cpp

test_core.exe: libcti_core.a ../encrypt/libtomlight.a $(TEST_OBJS) tests/test_core.o
	$(CXX) -o $@ $^ $(CLIBS)

tests/testExchanger.o: libcti_core.a $(SS_HEADERS) tests/testExchanger.cpp

tests/test_exchanger.o: ../catch1/catch.hpp tests/test_exchanger.cpp

test_exchanger.exe: libcti_core.a ../encrypt/libtomlight.a tests/testExchanger.o tests/test_exchanger.o
	$(CXX) -o $@ $^ $(CLIBS)

#-----------------------------------------------------------------
# modal analysis tools
#----------------------------------------------------------------

modal.o:	modal.cpp MiscUtils.hpp ByteSwap.hpp
	$(CXX) $(CXXFLAGS) -I../../ping -c $< -o $@

modal.exe:	modal.o libcti_core.a ../encrypt/libtomlight.a
	$(CXX) -o $@ $^ $(CLIBS)  $(LAPACK_LIB) $(BLAS_LIB)

modal:
	make modal.exe

#-----------------------------------------------------------------
# imaging stuff
#----------------------------------------------------------------

CtiCanvas_s.o:	CTI.hpp CvImageMap.hpp IntFlag.hpp CtiCanvas.hpp CtiCanvas.cpp ImageMetadata.hpp PngDataChunk.hpp WriteImageData.hpp PlaneData.hpp

CtiCanvas.o:	CTI.hpp CvImageMap.hpp IntFlag.hpp CtiCanvas.hpp CtiCanvas.cpp ImageMetadata.hpp PngDataChunk.hpp WriteImageData.hpp PlaneData.hpp

LesImageMapper_s.o: LesImageMapper.hpp LesImageMapper.cpp

LesImageMapper.o: LesImageMapper.hpp LesImageMapper.cpp

ColorMap_s.o: ColorMap.hpp ColorMap.cpp

ColorMap.o: ColorMap.hpp ColorMap.cpp

#-----------------------------------------------------------------
# solver stuff
#----------------------------------------------------------------

StaticSolver_s.o: $(SS_HEADERS) MiscUtils_s.o CtiRegister_s.o WebUI_s.o StaticSolver.cpp

StaticSolver.o: $(SS_HEADERS) MiscUtils.o CtiRegister.o WebUI.o StaticSolver.cpp LpData.hpp

StaticSolver_probes_s.o: $(SS_HEADERS) MiscUtils_s.o StaticSolver_probes.cpp

StaticSolver_probes.o: $(SS_HEADERS) MiscUtils.o StaticSolver_probes.cpp

StaticSolver_write_data_s.o: $(SS_HEADERS) MiscUtils_s.o StaticSolver_write_data.cpp

StaticSolver_write_data.o: $(SS_HEADERS) MiscUtils.o StaticSolver_write_data.cpp

StaticSolver_multigrid_s.o: $(SS_HEADERS) MiscUtils_s.o Multigrid_s.o StaticSolver_multigrid.cpp

StaticSolver_multigrid.o: $(SS_HEADERS) MiscUtils.o Multigrid.o StaticSolver_multigrid.cpp

# --------------------------------------------

../encrypt/libtomlight.a:
	make -C ../encrypt

# --------------------------------------------

clean:
	rm -f *.o *.a tests/*.o \
	testStripedMesh.exe testStaticSolver.exe testRegister.exe testMesh.exe testCti.exe test_core.exe

.PHONY: test testStripedMesh testStaticSolver testRegister modal testMesh testCti testSurfaceShm
