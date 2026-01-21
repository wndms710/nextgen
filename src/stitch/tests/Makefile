SHELL := /bin/bash
include ../../../Makefile.in

MPIRUN 	?=	mpirun -np 2
VALGRIND ?=	valgrind

include ../../../Makefile.in

# includes for the core (needs to be updated so its the include dir) 
CXXFLAGS +=	-I../../core/common
CXXFLAGS +=	-I../../core/encrypt

STITCH_OPT_FLAGS ?=	
CXXFLAGS +=	$(STITCH_OPT_FLAGS)

CTI_LIBS = 	-L../../core/common -lcti_core -L../../core/encrypt -ltomlight $(CLIBS)

default:	../stitch.exe ../stitch_s.exe zero.exe one.exe test1 test1x test2 test2x test3p test3s \
  		test_stitch_vg test_stitch_s_vg

tests_s:	../stitch_s.exe zero.exe one.exe test1 test1x test2 test2x test3s

../stitch.exe:
	make -C .. stitch

../stitch_s.exe:
	make -C .. stitch_s

zero.exe:
	$(CXX_NO_MPI) -o zero.exe ./zero.cpp

one.exe:
	$(CXX_NO_MPI) -o one.exe ./one.cpp

# here's an example of running an executable that returns 0...
test1:	zero.exe
	@/bin/echo -n "running test1 [zero]............................."
	@if [ $(shell ./zero.exe 2>&1 > log_zero && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_zero]"; fi

test1x:	zero.exe
	@/bin/echo -n "running test1x [zero checking for 1]............."
	@if [ $(shell ./zero.exe 2>&1 > log_zero && echo 0 || echo 1) = 1 ]; then echo "FAILED [see log_zero]"; else echo "PASSED"; fi

# here's an example of running an executable that returns 1...
test2:	one.exe
	@/bin/echo -n "running test2 [one].............................."
	@if [ $(shell ./one.exe 2>&1 > log_one && echo 0 || echo 1) = 1 ]; then echo "PASSED"; else echo "FAILED [see log_one]"; fi

test2x:	one.exe
	@/bin/echo -n "running test2x [one checking for 0].............."
	@if [ $(shell ./one.exe 2>&1 > log_one && echo 0 || echo 1) = 0 ]; then echo "FAILED [see log_one]"; else echo "PASSED"; fi

# here is an example of diff'ing the output log...
test3p_run:	../stitch.exe
	@/bin/echo -n "running test3p [stitch.box.in]..................."
	@../stitch.exe -i stitch.box.in 2>&1 > log_stitch.box.in

test3p:	test3p_run
	@if [ $(shell diff -I "User:" -I "Organization:" -I "Serial Number:" -I "version" -I "CTI_LICENSE" -I "WebUI run log" -I "elapsed time" -I "setup time" -I "bbox time" -I "surface ray intersection time" -I "bloated tri time" -I "time since last check" -I "timestamp (GMT)" log_stitch.box.in log_stitch.box.in.ref 2>&1 > log_stitch.box.in.diff && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_stitch.box.in.diff]"; fi

# here is a better example using serial stitch: you need to separate the 
# run from the diff utherwise make will evaluate the diff BEFORE the stitch_s.exe
# output is colected...
test3s_run:	../stitch_s.exe
	@/bin/echo -n "running test3s [stitch.box.in]..................."
	@../stitch_s.exe -i stitch.box.in 2>&1 > log_stitch.box.in_s

test3s:	test3s_run
	@if [ $(shell diff -I "User:" -I "Organization:" -I "Serial Number:" -I "version" -I "CTI_LICENSE" -I "WebUI run log" -I "elapsed time" -I "setup time" -I "bbox time" -I "surface ray intersection time" -I "bloated tri time" -I "time since last check" -I "timestamp (GMT)" log_stitch.box.in_s log_stitch.box.in.ref_s 2>&1 > log_stitch.box.in.diff_s && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_stitch.box.in.diff_s]"; fi

test_stitch_s_vg_run: ../stitch_s.exe
	@/bin/echo -n "running test_stitch_s_vg [test_stitch_vg.in]....."
	@$(VALGRIND) ../stitch_s.exe -i test_stitch_vg.in > log_test_stitch_s_vg 2>&1 

test_stitch_s_vg:	test_stitch_s_vg_run
	@if [ $(shell diff <(grep ' lost: ' log_test_stitch_s_vg | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_stitch_s_vg && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_stitch_s_vg]"; fi

test_stitch_vg_run: ../stitch.exe
	@/bin/echo -n "running test_stitch_vg [test_stitch_vg.in]......."
	@$(MPIRUN) $(VALGRIND) ../stitch.exe -i test_stitch_vg.in > log_test_stitch_vg 2>&1 

test_stitch_vg:	test_stitch_vg_run
	@if [ $(shell diff <(grep ' lost: ' log_test_stitch_vg | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_stitch_vg && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_stitch_vg]"; fi

# ==================================------===============================

clean:
	rm -vf *.o 
	rm -vf one.exe zero.exe
	rm -rf dump

.PHONY: test1 test1x test2 test2x test3 test3s test3s_run \
  test_stitch_vg test_stitch_vg_run test_stitch_s_vg test_stitch_s_vg_run
