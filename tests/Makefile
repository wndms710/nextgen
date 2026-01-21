

RUN=	mpirun -np 2 valgrind --leak-check=full --track-origins=yes

default:	build test_stitch test_charles test_helmholtz

build:
	@echo " Building executables"
	make -C ../src/surfer
	make -C ../src/stitch
	make -C ../src/charles charles

test_stitch:
	$(RUN) ../bin/stitch.exe -i stitch_test.in 2>&1 | tee log.stitch

test_charles:
	$(RUN) ../bin/charles.exe -i charles_test.in 2>&1 | tee log.charles

test_helmholtz:
	$(RUN) ../bin/charles.exe -i helmholtz_test.in 2>&1 | tee log.helmholtz


# catch1.x style tests
include ../Makefile.in

testPeriodicSector.o: testPeriodicSector.cpp

TEST_OBJECTS = testPeriodicSector.o \


test_runner.o: ../src/core/catch1/catch.hpp test_runner.cpp

test_runner.exe: $(TEST_OBJECTS) test_runner.o
	$(CXX) -o $@ $^


test_runner.exe:

clean:
	rm -rf dump
	rm -f log.stitch
	rm -f log.charles
	rm -f log.helmholtz
	rm -f result.sles
	rm -f *.o test_runner.exe

.PHONY:	build test_stitch test_charles test_helmholtz
