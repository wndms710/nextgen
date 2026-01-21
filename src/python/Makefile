# note: this is a local Makefile.in for now...
include ./Makefile.in
CXXFLAGS_PY += -DNO_MPI -I../chemistry -I../core/common -I../core/encrypt

default:
	@echo
	@echo "-------------------------------------------------"
	@echo "building cti python module..."
	@echo "-------------------------------------------------"
	@echo
	make build/lib.$(ARCH_PY)/cti.$(LIB_SUFFIX_PY)
	@echo
	@echo "-------------------------------------------------"
	@echo "Success! You can test the module as follows:"
	@echo "   cd build/lib.$(ARCH_PY)"
	@echo "   python (make sure you use the right path if you have multiple versions)"
	@echo "   >>> import cti"
	@echo "   >>> cti.chemtable_load(\"mychemtable\")"
	@echo
	@echo "To install:"
	@echo "   make install"
	@echo "-------------------------------------------------"
	@echo

PY_OBJS = ../core/common/NoMpi_py.o \
	../core/common/MpiStuff_py.o \
	../core/common/Guardian_py.o \
	../core/common/MiscUtils_py.o \
	../core/common/Params_py.o \
	../core/common/CtiRegister_py.o \
	../core/common/CTI_py.o \
	../core/encrypt/tomcrypt_py.o \
	../core/encrypt/tommath_py.o \
	../core/common/Version_py.o \
	../core/common/RestartHashUtilities_py.o \
	../core/common/Logger_py.o \
	ctimodule_py.o

%_py.o: %.cpp
	$(CXX_PY) $(CXXFLAGS_PY) -c $< -o $@

build/lib.$(ARCH_PY)/cti.$(LIB_SUFFIX_PY):	$(PY_OBJS)
	mkdir -p build/lib.$(ARCH_PY)
	$(CXX_PY) $^ $(LDFLAGS_PY) -o $@

# -------------------------------------------------------------

install:	build/lib.$(ARCH_PY)/cti.$(LIB_SUFFIX_PY)
	$(PYTHON) setup.py install
	@echo
	@echo "-------------------------------------------------"
	@echo "Success!"
	@echo "-------------------------------------------------"
	@echo

# -------------------------------------------------------------

clean:
	rm -vf $(PY_OBJS)
	rm -fR build







#---------------------------------
# main solver driver
#---------------------------------
#
#CTIMODULE_OBJECTS_S = 		ctimodule_s.o $(CORE_DEPS_S)
#
#ctimodule_s.o:	ctimodule.cpp
#	$(CXX_NO_MPI) -DWITHOUT_PYTHON -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@
#	@#@echo "     CXX_S $<"
#
#ctimodule_s.exe:		$(CTIMODULE_OBJECTS_S)
#	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS_S)

