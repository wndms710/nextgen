SHELL := /bin/bash
include ../../../Makefile.in

default: fwh.exe test_fwh 

fwh.exe: ../fwh.exe
	cp ../fwh.exe ./.

../fwh.exe:
	make -C ..

test_fwh_run: ../fwh.exe
	@/bin/echo -n "running test_fwh [completion]..................."
	@(mpirun -np 2 valgrind ../fwh.exe -i test_fwh.in &> log_test_fwh)  && echo "PASSED" || echo "FAILED [see log_test_fwh]"
	@/bin/echo -n "running test_fwh [vg]..................."

test_fwh:	test_fwh_run
	@if [ $(shell diff <(grep ' lost: ' log_test_fwh | awk '{ sum += $$4 } END { print sum }') <(echo "0") &>> log_test_fwh && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_fwh]"; fi


clean:
	rm -fv ./log*
	rm -fv ./*.mles
	rm -fv ./*.sles
	rm -fv ./*.exe
	rm -rfv ./dump

.PHONY: test_fwh 
