include ../../../Makefile.in

CXXFLAGS 	+= 	-I../../core/common
CXXFLAGS 	+= 	-I../../core/encrypt
CXXFLAGS 	+=	-I../

CXXFLAGS_NO_MPI +=	-I../../core/common
CXXFLAGS_NO_MPI +=	-I../../core/encrypt
CXXFLAGS_NO_MPI +=	-I../


SERIAL_LIBS	=	../../core/common/libcti_core_s.a \
			../../core/encrypt/libtomlight_s.a


tests :		test_query

../../core/common/libcti_core_s.a:
	make -C ../../core/common libcti_core_s

../../core/encrypt/libtomlight_s.a:
	make -C ../../core/encrypt libtomlight_s

%_s.o : %.cpp
	$(Q)$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI $< -c -o $@
	@echo "     CXX_S $<"

%_s.exe : $(SERIAL_LIBS) %_s.o
	$(CXX_NO_MPI) -o $@ $^ $(SERIAL_LIBS) $(CLIBS_NO_MPI) $(CANTERA_LIBS)



# --------------------------------
# simple test query of a chemtable
# --------------------------------

test_query:	test_query_s.exe
	@/bin/echo -n "running 2d chemtable query [completion]............................."
	@if [ $(shell ./test_query_s.exe >& log_query2d && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_zero]"; fi


