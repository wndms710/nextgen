SHELL := /bin/bash
include ../../../Makefile.in

MPIRUN 	?=	mpirun -np 2
VALGRIND ?=	valgrind

default: stitch.exe charles.exe test_runCheck test_lspDrag test_lspDrag_ell test_lpSolidDrag_1w test_lpSolidDrag_2w \
		test_lpLiquidDrag_1w test_lpLiquidDrag_2w test_charles_s_vg test_charles_vg \
  	test_helmholtz_s_vg test_helmholtz_vg test_alg_wm test_charles_fwh

stitch.exe: ../../stitch/stitch.exe
	cp ../../stitch/stitch.exe ./. 
	
charles.exe: ../charles.exe
	cp ../charles.exe ./.

../../stitch/stitch.exe:
	make -C ../../stitch

../charles.exe:
	make -C ..

../charles_s.exe:
	make -C .. charles_s

test_runCheck: charles.exe box.mles
	@/bin/echo -n "running test_runCheck [check successful run].............."
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lsp_runCheck.in 2>&1 > log_lsp_runCheck) && echo "PASSED" || echo "FAILED see log_lsp_runCheck"

box.mles: stitch.exe 
	@/bin/echo -n "making an HCP grid........................................"
	@($(MPIRUN) ./stitch.exe -i ./assets/stitch.box.in 2>&1 > log_stitch_box) && echo "PASSED" || "FAILED see log_stitch_box"

test_lspDrag: charles.exe CompData_lsp_drag.exe box.mles 
	@/bin/echo -n "running test_lspDrag [test particle deceleration]........."
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lsp_drag.in 2>&1 > log_lsp_drag) && (grep 'starting step' log_lsp_drag | sed  's/.*time://' | sed 's/dt:.*//' > file1) && (grep 'dumpRange: up,' log_lsp_drag | sed 's/.*0://' | sed 's/:.*//1' > file2) && (paste -d " " file1 file2 > log_lsp_drag_vel.dat) && (rm file1 file2) && ./CompData_lsp_drag.exe -i log_lsp_drag_vel.dat -r ./assets/reference_lsp_drag.dat -t 1e-2 2>&1 >> log_lsp_drag_vel.dat && echo "PASSED" || echo "FAILED see log_lsp_drag and log_lsp_drag_vel.dat"

test_lspDrag_ell: charles.exe CompData_lsp_drag.exe box.mles 
	@/bin/echo -n "running test_lspDrag_ell [test ellipsoidal particle deceleration] "
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lsp_drag_ell.in 2>&1 > log_lsp_drag_ell) && (grep 'starting step' log_lsp_drag_ell | sed  's/.*time://' | sed 's/dt:.*//' > file1) && (grep 'dumpRange: up,' log_lsp_drag_ell | sed 's/.*0://' | sed 's/:.*//1' > file2) && (paste -d " " file1 file2 > log_lsp_drag_vel_ell.dat) && (rm file1 file2) && ./CompData_lsp_drag.exe -i log_lsp_drag_vel_ell.dat -r ./assets/reference_lsp_drag.dat -t 1e-2 2>&1 >> log_lsp_drag_vel_ell.dat && echo "PASSED" || echo "FAILED see log_lsp_drag_ell and log_lsp_drag_vel_ell.dat"

test_lpSolidDrag_1w: charles.exe CompData_lsp_drag.exe box.mles 
	@/bin/echo -n "running test_lpSolidDrag_1w [test one-way lp Solid]......."
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lpSolid_drag_1w.in 2>&1 > log_lpSolid_drag_1w) && (grep 'starting step' log_lpSolid_drag_1w | sed  's/.*time://' | sed 's/dt:.*//' > file1) && (grep 'dumpRange: up,' log_lpSolid_drag_1w | sed 's/.*0://' | sed 's/:.*//1' > file2) && (paste -d " " file1 file2 > log_lpSolid_drag_vel_1w.dat) && (rm file1 file2) && ./CompData_lsp_drag.exe -i log_lpSolid_drag_vel_1w.dat -r ./assets/reference_lpSolid_drag.dat -t 1e-2 2>&1 >> log_lpSolid_drag_vel_1w.dat && echo "PASSED" || echo "FAILED see log_lpSolid_drag_1w and log_lpSolid_drag_vel_1w.dat"

test_lpSolidDrag_2w: charles.exe CompData_lsp_drag.exe box.mles 
	@/bin/echo -n "running test_lpSolidDrag_2w [test two-way lp solid]......."
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lpSolid_drag_2w.in 2>&1 > log_lpSolid_drag_2w) && (grep 'starting step' log_lpSolid_drag_2w | sed  's/.*time://' | sed 's/dt:.*//' > file1) && (grep 'dumpRange: up,' log_lpSolid_drag_2w | sed 's/.*0://' | sed 's/:.*//1' > file2) && (paste -d " " file1 file2 > log_lpSolid_drag_vel_2w.dat) && (rm file1 file2) && ./CompData_lsp_drag.exe -i log_lpSolid_drag_vel_2w.dat -r ./assets/reference_lpSolid_drag.dat -t 1e-2 2>&1 >> log_lpSolid_drag_vel_2w.dat && echo "PASSED" || echo "FAILED see log_lpSolid_drag_2w and log_lpSolid_drag_vel_2w.dat"

test_lpLiquidDrag_1w: charles.exe CompData_lsp_drag.exe box.mles 
	@/bin/echo -n "running test_lpLiquidDrag_1w [test one-way lp liquid]....."
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lpLiquid_drag_1w.in 2>&1 > log_lpLiquid_drag_1w) && (grep 'starting step' log_lpLiquid_drag_1w | sed  's/.*time://' | sed 's/dt:.*//' > file1) && (grep 'dumpRange: up,' log_lpLiquid_drag_1w | sed 's/.*0://' | sed 's/:.*//1' > file2) && (paste -d " " file1 file2 > log_lpLiquid_drag_vel_1w.dat) && (rm file1 file2) && ./CompData_lsp_drag.exe -i log_lpLiquid_drag_vel_1w.dat -r ./assets/reference_lpLiquid_drag.dat -t 1e-2 2>&1 >> log_lpLiquid_drag_vel_1w.dat && echo "PASSED" || echo "FAILED see log_lpLiquid_drag_1w and log_lpLiquid_drag_vel_1w.dat"

test_lpLiquidDrag_2w: charles.exe CompData_lsp_drag.exe box.mles 
	@/bin/echo -n "running test_lpLiquidDrag_2w [test two-way lp liquid]....."
	@($(MPIRUN) ./charles.exe -i ./assets/charles_lpLiquid_drag_2w.in 2>&1 > log_lpLiquid_drag_2w) && (grep 'starting step' log_lpLiquid_drag_2w | sed  's/.*time://' | sed 's/dt:.*//' > file1) && (grep 'dumpRange: up,' log_lpLiquid_drag_2w | sed 's/.*0://' | sed 's/:.*//1' > file2) && (paste -d " " file1 file2 > log_lpLiquid_drag_vel_2w.dat) && (rm file1 file2) && ./CompData_lsp_drag.exe -i log_lpLiquid_drag_vel_2w.dat -r ./assets/reference_lpLiquid_drag.dat -t 1e-2 2>&1 >> log_lpLiquid_drag_vel_2w.dat && echo "PASSED" || echo "FAILED see log_lpLiquid_drag_2w and log_lpLiquid_drag_vel_2w.dat"

CompData_lsp_drag.exe: CompData.o CompData_lsp_drag.cpp
	$(CXX_NO_MPI) CompData_lsp_drag.cpp -o $@ $<

CompData.o: CompData.cpp CompData.hpp
	$(CXX_NO_MPI) -c $<

test_charles_s_vg_run: ../charles_s.exe
	@/bin/echo -n "running test_charles_s_vg [test_charles_vg.in]............"
	@$(VALGRIND) ../charles_s.exe -i test_charles_vg.in > log_test_charles_s_vg 2>&1 

test_charles_s_vg:	test_charles_s_vg_run
	@if [ $(shell diff <(grep ' lost: ' log_test_charles_s_vg | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_charles_s_vg && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_charles_s_vg]"; fi

test_charles_vg_run: ../charles.exe
	@/bin/echo -n "running test_charles_vg [test_charles_vg.in].............."
	@$(MPIRUN) $(VALGRIND) ../charles.exe -i test_charles_vg.in > log_test_charles_vg  2>&1 

test_charles_vg:	test_charles_vg_run
	@if [ $(shell diff <(grep ' lost: ' log_test_charles_vg | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_charles_vg && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_charles_vg]"; fi

test_helmholtz_s_vg_run: ../charles_s.exe
	@/bin/echo -n "running test_helmholtz_s_vg [test_helmholtz_vg.in]........"
	@$(VALGRIND) ../charles_s.exe -i test_helmholtz_vg.in > log_test_helmholtz_s_vg 2>&1 

test_helmholtz_s_vg:	test_helmholtz_s_vg_run
	@if [ $(shell diff <(grep ' lost: ' log_test_helmholtz_s_vg | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_helmholtz_s_vg && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_helmholtz_s_vg]"; fi

test_helmholtz_vg_run: ../charles.exe
	@/bin/echo -n "running test_helmholtz_vg [test_helmholtz_vg.in].........."
	@$(MPIRUN) $(VALGRIND) ../charles.exe -i test_helmholtz_vg.in > log_test_helmholtz_vg 2>&1 

test_helmholtz_vg:	test_helmholtz_vg_run
	@if [ $(shell diff <(grep ' lost: ' log_test_helmholtz_vg | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_helmholtz_vg && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_helmholtz_vg]"; fi

wm_test.exe:	wm_test.cpp
	$(CXX_NO_MPI) -I../wm -o $@ $<

test_alg_wm:	wm_test.exe
	@/bin/echo -n "running algebraic wm test ................................"
	@if [ $(shell ./wm_test.exe >& log_wm_test && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_wm_test]"; fi

test_charles_fwh_run: ../charles.exe
	@/bin/echo -n "running test_charles_fwh [completion]....................."
	@($(MPIRUN) $(VALGRIND) ../charles.exe -i test_charles_fwh.in > log_test_charles_fwh 2>&1)  && echo "PASSED" || echo "FAILED [see log_test_charles_fwh]"
	@/bin/echo -n "running test_charles_fwh [vg]............................."

test_charles_fwh:	test_charles_fwh_run
	@if [ $(shell diff <(grep ' lost: ' log_test_charles_fwh | awk '{ sum += $$4 } END { print sum }') <(echo "0") 2>&1 >> log_test_charles_fwh && echo 0 || echo 1) = 0 ]; then echo "PASSED"; else echo "FAILED [see log_test_charles_fwh]"; fi


clean:
	rm -fv ./log*
	rm -fv ./*.mles
	rm -fv ./*.sles
	rm -fv ./*.exe
	rm -rfv ./dump
	rm -fv ./*.o

.PHONY: test_runCheck test_lspDrag test_lpSolidDrag_1w test_lpSolidDrag_2w test_lpLiquidDrag_1w test_lpLiquidDrag_2w test_lspDrag_ell\
	test_charles_s_vg_run test_charles_s_vg test_charles_vg_run test_charles_vg \
	test_helmholtz_s_vg_run test_helmholtz_s_vg test_helmholtz_vg_run test_helmholtz_vg \
	test_alg_wm
