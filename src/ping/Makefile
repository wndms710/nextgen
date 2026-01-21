include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir)
CXXFLAGS_NO_MPI +=	-I../core/common
CXXFLAGS_NO_MPI +=	-I../core/encrypt

CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt

CXXFLAGS +=	-O3

CTI_LIBS = 	-L../core/common -lcti_core_s -L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

# optional flags and parameters to set if you are using
# imageDmd (which still exists in serial...)
CXXFLAGS        +=     -DCACHE_IMAGES
CXXFLAGS_NO_MPI +=     -DCACHE_IMAGES

#LAPACK_LIB ?=		-llapack
#BLAS_LIB ?=		-lblas

ifndef LAPACK_LIB
default: ping
else ifndef BLAS_LIB
default: ping
else
default: ping modal
endif


../../bin/%.exe : %.exe
	cp $^ ../../bin/

ping:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for ping..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making ping_s.exe..."
	@echo "-------------------------------------------------"
	@echo
	make ping_s.exe
	make ../../bin/ping_s.exe

modal:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for modal suite..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt libtomlight.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making imagePod.exe..."
	@echo "-------------------------------------------------"
	@echo
	make imagePod.exe
	make ../../bin/imagePod.exe
	@echo
	@echo "-------------------------------------------------"
	@echo "making imageDmd.exe..."
	@echo "-------------------------------------------------"
	@echo
	make imageDmd.exe
	make ../../bin/imageDmd.exe

imageDiff:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for image diff utility..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making imageDiff.exe..."
	@echo "-------------------------------------------------"
	@echo
	make imageDiff.exe

imageDmd:
	make imageDmd.exe

imagePod:
	make imagePod.exe



PNG_DEPS = PngData.hpp PngCache.hpp ImageTools.hpp CtiType.hpp Journal.hpp ImageLegend.hpp
POD_DEPS = ../core/common/Lapack.hpp Pod.hpp PngAnalysis.hpp PngData.hpp ImageLegend.hpp
DMD_DEPS = ../core/common/Lapack.hpp Dmd.hpp PngAnalysis.hpp PngData.hpp ImageLegend.hpp
# PngImage.hpp, ImageMetaData.hpp PngDataChunk.hpp, ImageMetaData.hpp and ColorMap.hpp are also dependencies, but are in ../core/common
Journal_s.o: 	Journal.hpp Journal.cpp

imagePod.o:	$(POD_DEPS) imagePod.cpp
imageDmd.o:	$(DMD_DEPS) imageDmd.cpp
imageDmdCompositeNew.o:	$(HDR_DEPS) imageDmdCompositeNew.cpp
imageProbe.o:	$(HDR_DEPS) imageProbe.cpp
rescale.o:	$(HDR_DEPS) rescale.cpp

ping_s.o:	Journal_s.o $(PNG_DEPS) ping.cpp
imagePod_s.o:	$(POD_DEPS) imagePod.cpp
imageDmd_s.o:	$(DMD_DEPS) imageDmd.cpp
imageProbe_s.o:	$(HDR_DEPS) imageProbe.cpp
rescale_s.o:	$(HDR_DEPS) rescale.cpp

../core/common/libcti_core.a : 
	make -C  ../core/common libcti_core.a

../core/common/libcti_core_s.a : 
	make -C  ../core/common libcti_core_s.a

ping_s.exe :	../core/common/libcti_core_s.a Journal_s.o ping_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)

imagePod.exe : imagePod.o ../core/common/libcti_core.a 
	$(CXX) -o $@ $^ $(LAPACK_LIB) $(BLAS_LIB) $(CLIBS) 

imagePod_s.exe: imagePod_s.o ../core/common/libcti_core_s.a 
	$(CXX_NO_MPI) -o $@ $^ $(LAPACK_LIB) $(BLAS_LIB) $(CLIBS) 

imageDmd.exe : imageDmd.o ../core/common/libcti_core.a 
	$(CXX) -o $@ $^ $(LAPACK_LIB) $(BLAS_LIB) $(CLIBS) 

imageDmd_s.exe : imageDmd_s.o ../core/common/libcti_core_s.a 
	$(CXX_NO_MPI) -o $@ $^ $(LAPACK_LIB) $(BLAS_LIB) $(CLIBS) 

imageDmdCompositeNew.exe : imageDmdCompositeNew.o
	$(CXX) -o $@ $^ $(LAPACK_LIB) $(BLAS_LIB) $(CLIBS)

imageDiff_s.o: $(PNG_DEPS) imageDiff.cpp

imageDiff.exe: ../core/common/libcti_core_s.a imageDiff_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)


%_s.o : %.cpp
	$(Q)$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI $< -c -o $@
	@echo "     CXX_S $<"

clean:
	rm -vf *.o
	rm -vf *.a
	rm -vf *.exe
	make -C test veryclean

.PHONY: imageTools imageDiff modal ping imageDmd imagePod
