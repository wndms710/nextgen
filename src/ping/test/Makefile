include ../../../Makefile.in

SHELL=/bin/bash

# includes for the core (needs to be updated so its the include dir)
CXXFLAGS_NO_MPI +=	-I../../core/common
CXXFLAGS_NO_MPI +=	-I../../core/encrypt

#########CTI_LIBS = 	-L../core/common -lcti_core_s  -L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

%_s.o : %.cpp
	$(Q)$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -DNO_MPI $< -c -o $@
	@echo "     CXX_S $<"

default:	tests_s

tests_s:	writeImageSingle_s writeImageTemporal_s probeImage_s

../ping_s.exe:
	make -C ../ ping

../../core/common/libcti_core_s.a:
	make -C ../../core/common libcti_core_s.a

stitch.exe: ../../stitch/stitch.exe
	cp ../../stitch/stitch.exe ./. 
	
charles.cart_canvas.exe: assets/charles.exe
	cp assets/charles.exe ./charles.cart_canvas.exe

../../stitch/stitch.exe:
	make -C ../../stitch

assets/charles.exe:
	make -C assets/ charles

cart_canvas.mles: stitch.exe 
	@/bin/echo -n "making a cart grid for imaging................................................................."
	@($(MPIRUN) ./stitch.exe -i ./assets/stitch.cart_canvas.in 2>&1 > log_stitch_cart_canvas) && echo "PASSED" || "FAILED see log_stitch_cart_canvas"

images_probe_test/phi.00000001.png: cart_canvas.mles charles.cart_canvas.exe                                               
	@/bin/echo -n "imaging a known variable on cart_canvas.mles..................................................."
	@($(MPIRUN) ./charles.cart_canvas.exe -i ./assets/charles.cart_canvas.in 2>&1 > log_charles_cart_canvas)
	@if [ -f images_probe_test/phi.00000001.png ]; then echo "PASSED"; else echo "FAILED [see log_charles_cart_canvas]"; fi

probeImage_s: images_probe_test/phi.00000001.png ../ping_s.exe
	@/bin/echo -n "tesing that probe of rgb value 255 (with RANGE 0:256) returns 254.996.........................."
	@../ping_s.exe --IMAGES images_probe_test/phi.00000001.png --PROBE GEOM POINT IJ 775 775 images_probe_test/probe_255.dat >> log_probeImage
	@val=$$(awk 'NR == 4 {print $$9}' images_probe_test/probe_255.dat); if [ "$$val" = "254.996" ]; then echo "PASSED"; else echo "FAILED val=$$val [see log_probeImage]"; fi

# tests
writeImage_s.exe: writeImage.cpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) -o writeImage_s.exe ./writeImage.cpp -L../../core/common -lcti_core_s $(CLIBS_NO_MPI)
writeImageSingle_s_run:	../../core/common/libcti_core_s.a writeImage_s.exe
			@/bin/echo -n "running writeImage --SINGLE...................................................................."
			@./writeImage_s.exe --SINGLE test.png >& log_writeImage_SINGLE 
writeImageSingle_s:	writeImageSingle_s_run
			@if [ -f test.png ]; then echo "PASSED"; else echo "FAILED [see log_writeImage_SINGLE]"; fi

writeImageTemporal_s_run:	../../core/common/libcti_core_s.a ../ping_s.exe writeImage_s.exe
			@/bin/echo -n "running writeImage --TEMPORAL.................................................................."
			@./writeImage_s.exe --TEMPORAL test.png >& log_writeImage_TEMPORAL  #generate reference images
			@../ping_s.exe --IMAGES images_ref/temporal_??.png --TEMPORAL_STATS >> log_writeImage_TEMPORAL 2>&1  #compute stats on images
			@../ping_s.exe --IMAGES images_ref/stats_avg.png images_stats/temporal_00.mean.png --DIFF_STATS >> log_writeImage_TEMPORAL 2>&1  #diff avg images
			@../ping_s.exe --IMAGES images_ref/stats_rms.png images_stats/temporal_00.rms.png --DIFF_STATS >> log_writeImage_TEMPORAL 2>&1  #diff avg images
			@../ping_s.exe --IMAGES images_ref/stats_min.png images_stats/temporal_00.min.png --DIFF_STATS >> log_writeImage_TEMPORAL 2>&1  #diff avg images
			@../ping_s.exe --IMAGES images_ref/stats_max.png images_stats/temporal_00.max.png --DIFF_STATS >> log_writeImage_TEMPORAL 2>&1  #diff avg images
writeImageTemporal_s:	writeImageTemporal_s_run #count diff statements in log to verify images match
			@count=$$(($$(grep "rgb 0(0)" log_writeImage_TEMPORAL | wc -l)+$$(grep "daTa 0(0)" log_writeImage_TEMPORAL | wc -l)+$$(grep "dpTh 0(0)" log_writeImage_TEMPORAL | wc -l)+$$(grep "zoNe 0(0)" log_writeImage_TEMPORAL | wc -l))); count_metadata=$$(grep "metadata-discrepancy" log_writeImage_TEMPORAL | wc -l); if [ "$$count" -eq "16" -a "$$count_metadata" -eq "0" ]; then echo "PASSED"; else echo "FAILED [see log_writeImage_TEMPORAL]"; fi


testCtiType_s.exe: testCtiType.cpp
	$(CXX_NO_MPI) $(CXXFLAGS_NO_MPI) -o testCtiType_s.exe ./testCtiType.cpp -L../../core/common -lcti_core_s $(CLIBS_NO_MPI)

# ==================================------===============================

clean:
	rm -vf *.o 
	rm -vf writeImage_s.exe 
	rm -vf testCtiType_s.exe
	rm -vf *.png
	rm -vf log_*
	rm -vrf images_ref
	rm -vrf images_stats
	rm -vf cart_canvas.mles
	rm -vrf images_probe_test

veryclean: clean
	rm -vf stitch.exe
	rm -vf charles.cart_canvas.exe
	rm -vf assets/charles.exe
	rm -vf assets/charles.o


.PHONY: tests_s writeImage_s writeImage_s_run probeImage_s


