
include ../../Makefile.in

default:       	lesinfo

CXXFLAGS_NO_MPI= 	-O2 -I../core/encrypt -I../core/common
CLIBS_NO_MPI+=		-L../core/common -lcti_core_s -L../core/encrypt/ -ltomlight_s 

../../bin/%.exe : %.exe
	cp $^ ../../bin/

lesinfo:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components lesinfo..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making lesinfo..."
	@echo "-------------------------------------------------"
	@echo
	make lesinfo.exe
	make ../../bin/lesinfo.exe


lesinfo.exe:		lesinfo.cpp 
	$(CXX_NO_MPI) -o $@ -DNO_MPI $(CXXFLAGS_NO_MPI) $^ $(CLIBS_NO_MPI)

lesdiff:	
	make lesdiff.exe
	make ../../bin/lesdiff.exe

lesdiff.exe:		lesdiff.cpp
	$(CXX_NO_MPI) -o $@ -DNO_MPI $(CXXFLAGS_NO_MPI) $^ $(CLIBS_NO_MPI)

# =================================================================

hold:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components hold..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common Logger.o
	@echo
	@echo "-------------------------------------------------"
	@echo "making hold..."
	@echo "-------------------------------------------------"
	@echo
	make hold.exe
	make ../../bin/hold.exe

hold.exe: 		hold.cpp 
	$(CXX) -o $@ $(CXXFLAGS) -I../core/common/ $^ ../core/common/Logger.o


# =================================================================

clean:
	rm -f *.exe core core.* 
	rm -f ../../bin/lesinfo.exe ../../bin/hold.exe
	rm -fR ii_files *~

.PHONY:	lesinfo hold lesdiff
