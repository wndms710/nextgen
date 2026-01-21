include ../../../Makefile.in

CXXFLAGS += 	-I../../core/common
CXXFLAGS +=	-I../../core/encrypt
CXXFLAGS += 	-I../../chemistry
CXXFLAGS += 	-I../ideal_gas
CXXFLAGS += 	-I../
CXXFLAGS += 	-I./

CHARLES_OPT_FLAGS ?=	
CXXFLAGS +=	$(CHARLES_OPT_FLAGS)

cti_lsp:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/charles/lsp/libcti_lsp.a..."
	@echo "-------------------------------------------------"
	@echo
	make libcti_lsp.a

CtiLiquid.o: CtiLiquid.hpp CtiLiquid.cpp LpStuff.hpp LpContainer.hpp Lsp.hpp
CtiSolid.o: CtiSolid.hpp CtiSolid.cpp LpStuff.hpp LpContainer.hpp Lsp.hpp
LpInjector.o: LpInjector.hpp LpInjector.cpp LpStuff.hpp LpContainer.hpp Lsp.hpp
LspInjector.o: LspInjector.hpp LspInjector.cpp LpStuff.hpp LpContainer.hpp Lsp.hpp
LspStats.o: LspStats.hpp LspStats.cpp LpStuff.hpp LpContainer.hpp Lsp.hpp
LspBcs.o: LspBcs.hpp LpStuff.hpp LpContainer.hpp Lsp.hpp

LSP_OBJS= CtiLiquid.o CtiSolid.o LspStats.o LpInjector.o LspInjector.o LspBcs.o
libcti_lsp.a: $(LSP_OBJS)
	$(AR) $@ $^
	$(RANLIB) $@

# ==================================------===============================
# serial lsp...

CXXFLAGS_NO_MPI += 	-I../../core/common
CXXFLAGS_NO_MPI +=	-I../../core/encrypt
CXXFLAGS_NO_MPI += 	-I../../chemistry
CXXFLAGS_NO_MPI	+= 	-I../ideal_gas
CXXFLAGS_NO_MPI += 	-I../
CXXFLAGS_NO_MPI += 	-I./

CXXFLAGS_NO_MPI +=	$(CHARLES_OPT_FLAGS)

#Rule for serial compilation
%_s.o : %.cpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@
	@#@echo "     CXX_S $<"

cti_lsp_s:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/charles/lsp/libcti_lsp_s.a..."
	@echo "-------------------------------------------------"
	@echo
	make libcti_lsp_s.a

CtiLiquid_s.o: CtiLiquid.hpp CtiLiquid.cpp LpStuff.hpp Lsp.hpp
CtiSolid_s.o: CtiSolid.hpp CtiSolid.cpp LpStuff.hpp Lsp.hpp
LpInjector_s.o: LpInjector.hpp LpInjector.cpp LpStuff.hpp Lsp.hpp
LspInjector_s.o: LspInjector.hpp LspInjector.cpp LpStuff.hpp Lsp.hpp
LspStats_s.o: LspStats.hpp LspStats.cpp LpStuff.hpp Lsp.hpp
LspBcs_s.o: LspBcs.hpp LpStuff.hpp Lsp.hpp

LSP_OBJS_S= CtiLiquid_s.o CtiSolid_s.o LspStats_s.o LpInjector_s.o LspInjector_s.o LspBcs_s.o
libcti_lsp_s.a: $(LSP_OBJS_S)
	$(AR) $@ $^
	$(RANLIB) $@

# ==================================------===============================

CTI_LIBS = -L. -lcti_lsp -L../../core/common -lcti_core -L../../core/encrypt -ltomlight
testCtiLiquid:
	@echo
	@echo "-------------------------------------------------"
	@echo "making src/charles/lsp/testCtiLiquid.exe..."
	@echo "-------------------------------------------------"
	@echo
	make libcti_lsp.a
	make testCtiLiquid.exe
testCtiLiquid.o: testCtiLiquid.cpp CtiLiquid.cpp CtiLiquid.hpp
testCtiLiquid.exe: testCtiLiquid.o
	$(CXX) -o $@ $^ $(CTI_LIBS)

#============================================

clean:
	rm -vf *.o *.exe
	rm -vf libcti_lsp.a
	rm -vf libcti_lsp_s.a
