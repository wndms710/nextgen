include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir) 
CXXFLAGS +=     -I../surfer2
CXXFLAGS +=     -I../stitch2
CXXFLAGS +=	-I../core/hdrlib
CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt

CTI_LIBS = 	-L../surfer2 -lcti_image2 \
		-L../core/common -lcti_core2 \
		-L../core/encrypt -ltomlight $(CLIBS)

#CTI_LIBS = 	../surfer2/CtiCanvas.o ../surfer2/LesImageMapper.o \
#		-L../core/common -lcti_core2 \
#		-L../core/encrypt -ltomlight $(CLIBS)

# includes for the core (needs to be updated so its the include dir)
CXXFLAGS_NO_MPI +=	-I../core/common
CXXFLAGS_NO_MPI +=	-I../core/encrypt

CTI_LIBS_NO_MPI = -L../core/common -lcti_core_s -L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

#default:	postpro
default:	lespost

postpro:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for postpro..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core2.a
	make -C ../surfer2 libcti_image2.a 
	make -C ../core/encrypt 
	@echo
	@echo "-------------------------------------------------"
	@echo "making postpro.exe..."
	@echo "-------------------------------------------------"
	@echo
	make postpro.exe

CORE_HDR_DEPS=	../core/common/StaticSolver.hpp ../core/common/StripedMesh.hpp

postpro.o :	../core/common/libcti_core2.a $(CORE_HDR_DEPS) PostproSolver.hpp postpro.cpp

postpro.exe :	postpro.o
	$(CXX) -o $@ $^ $(CTI_LIBS)

#-------------------------------------

lespost:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for lespost..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making lespost.exe..."
	@echo "-------------------------------------------------"
	@echo
	make lespost.exe
#	make ../../bin/surfer.exe

lespost_s.o:	lespost.cpp ../core/common/libcti_core_s.a LesPost.hpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@

lespost.exe:	lespost_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS_NO_MPI)


#-------------------------------------

ws_post.exe: 	ws_post.cpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -o $@

#-------------------------------------

#Rule for serial compilation
#%_s.o : %.cpp
#	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@
#	@#@echo "     CXX_S $<"

#-------------------------------------

clean:
	rm -vf *.o 
	rm -vf postpro.exe lespost.exe ws_post.exe

.PHONY: postpro lespost
