include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir)
CXXFLAGS 		+= -I../core/common
CXXFLAGS 		+= -I../core/encrypt
CXXFLAGS 		+= -I../chemistry
CXXFLAGS 		+= -I../quiver/charles
CXXFLAGS 		+= -I./lsp
CXXFLAGS 		+= -I./ideal_gas
CXXFLAGS 		+= -I./premixed
CXXFLAGS 		+= -I./nonpremixed
CXXFLAGS 		+= -I./helmholtz
CXXFLAGS 		+= -I.

CHARLES_OPT_FLAGS 	?=
CXXFLAGS 	  	+=	$(CHARLES_OPT_FLAGS)

CTI_LIBS = 		-L./lsp -lcti_lsp -L../core/common -lcti_core \
			-L../core/encrypt -ltomlight $(CLIBS)

default:	charles

../../bin/%.exe : %.exe
	cp $^ ../../bin/

charles:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles.exe
	make ../../bin/charles.exe


CORE_DEPS=			../core/common/libcti_core.a

FLOW_HDR_DEPS=			FlowSolver.hpp Cht.hpp lsp/LpStuff.hpp lsp/LpContainer.hpp lsp/LpInjector.hpp lsp/LpBcs.hpp lsp/CtiLiquid.hpp \
										lsp/CtiSolid.hpp lsp/CtiMaterial.hpp ../core/common/LpData.hpp

LSP_HDR_DEPS=			lsp/LpStuff.hpp lsp/LpContainer.hpp lsp/Lsp.hpp lsp/LpInjector.hpp lsp/LspInjector.hpp lsp/LspStats.hpp \
				lsp/CtiMaterial.hpp lsp/CtiLiquid.hpp lsp/CtiSolid.hpp lsp/LpInjector.cpp lsp/LspInjector.cpp lsp/LspStats.cpp \
				lsp/CtiLiquid.cpp lsp/CtiSolid.cpp lsp/LspBcs.cpp lsp/LspBcs.hpp lsp/LpBcs.hpp ../core/common/LpData.hpp

#LSP_SOLVER_DEPS=		lsp/IdealGasSolverWithLsp.hpp lsp/SingleLspSolver.hpp lsp/libcti_lsp.a
LSP_SOLVER_DEPS=		lsp/IdealGasSolverWithLsp.hpp lsp/libcti_lsp.a

lsp/libcti_lsp.a:       	$(LSP_HDR_DEPS)
	make -C ./lsp cti_lsp

#-----------------------------------------
# ideal gas solver
#-----------------------------------------

IG_OBJECTS=				ideal_gas/IdealGasSolver.o ideal_gas/IdealGasSolverSgs.o \
					ideal_gas/IdealGasSolverIo.o ideal_gas/IdealGasSolverBcs.o

ideal_gas/IdealGasSolver.o:		ideal_gas/IdealGasSolver.hpp ideal_gas/IdealGasSolverFlux.hpp \
					$(CORE_DEPS) $(FLOW_HDR_DEPS) ideal_gas/IdealGasSolver.cpp

ideal_gas/IdealGasSolverSgs.o:		ideal_gas/IdealGasSolver.hpp Sgs.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ideal_gas/IdealGasSolverSgs.cpp

ideal_gas/IdealGasSolverIo.o:		ideal_gas/IdealGasSolver.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ideal_gas/IdealGasSolverIo.cpp bcprofile/ProfileReader.hpp

ideal_gas/IdealGasSolverBcs.o:		ideal_gas/IdealGasSolverBcs.hpp ideal_gas/IdealGasSolver.hpp \
					ideal_gas/IdealGasSolverFlux.hpp \
					it/InflowTurbulenceIG.hpp bcprofile/ProfileReader.hpp BcCommon.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ideal_gas/IdealGasSolverBcs.cpp

charles_ig:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles_ig..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles_ig.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles_ig.exe
	make ../../bin/charles_ig.exe

CHARLES_IG_OBJECTS = 		charles_ig.o $(IG_OBJECTS)

charles_ig.o :			$(CORE_DEPS) $(FLOW_HDR_DEPS) $(LSP_SOLVER_DEPS) \
				States.hpp BcCommon.hpp CtiMemory.hpp CommContainers.hpp \
				ideal_gas/IdealGasSolver.hpp charles_ig.cpp

charles_ig.exe :		$(CHARLES_IG_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

#-----------------------------------------
# premixed solver
#-----------------------------------------

PREMIXED_OBJECTS=			premixed/PremixedSolver.o premixed/PremixedSolverSgs.o \
					premixed/PremixedSolverIo.o premixed/PremixedSolverBcs.o

premixed/PremixedSolver.o:		premixed/PremixedSolver.hpp premixed/PremixedSolverFlux.hpp \
					$(CORE_DEPS) $(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp  premixed/PremixedSolver.cpp

premixed/PremixedSolverSgs.o:		premixed/PremixedSolver.hpp Sgs.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp premixed/PremixedSolverSgs.cpp

premixed/PremixedSolverIo.o:		premixed/PremixedSolver.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp premixed/PremixedSolverIo.cpp

premixed/PremixedSolverBcs.o:		premixed/PremixedSolverBcs.hpp premixed/PremixedSolver.hpp \
					premixed/PremixedSolverFlux.hpp \
					it/InflowTurbulenceP.hpp bcprofile/ProfileReader.hpp BcCommon.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp premixed/PremixedSolverBcs.cpp

charles_pm:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles_pm..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles_pm.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles_pm.exe
	make ../../bin/charles_pm.exe

CHARLES_PM_OBJECTS = 		charles_pm.o $(PREMIXED_OBJECTS)

charles_pm.o :			$(CORE_DEPS) $(FLOW_HDR_DEPS) $(LSP_SOLVER_DEPS) \
				States.hpp BcCommon.hpp CtiMemory.hpp CommContainers.hpp \
				premixed/PremixedSolver.hpp charles_pm.cpp

charles_pm.exe :		$(CHARLES_PM_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

#-----------------------------------------
# nonpremixed solver
#-----------------------------------------

NONPREMIXED_OBJECTS=			nonpremixed/NonpremixedSolver.o nonpremixed/NonpremixedSolverSgs.o \
					nonpremixed/NonpremixedSolverIo.o nonpremixed/NonpremixedSolverBcs.o

nonpremixed/NonpremixedSolver.o:	nonpremixed/NonpremixedSolver.hpp premixed/PremixedSolverFlux.hpp \
					$(CORE_DEPS) $(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolver.cpp

nonpremixed/NonpremixedSolverSgs.o:	nonpremixed/NonpremixedSolver.hpp Sgs.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolverSgs.cpp

nonpremixed/NonpremixedSolverIo.o:	nonpremixed/NonpremixedSolver.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolverIo.cpp

nonpremixed/NonpremixedSolverBcs.o:	nonpremixed/NonpremixedSolverBcs.hpp nonpremixed/NonpremixedSolver.hpp \
					premixed/PremixedSolverFlux.hpp bcprofile/ProfileReader.hpp BcCommon.hpp $(CORE_DEPS) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolverBcs.cpp

charles_npm:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles_npm..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles_npm.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles_npm.exe
	make ../../bin/charles_npm.exe

CHARLES_NPM_OBJECTS = 		charles_npm.o $(NONPREMIXED_OBJECTS)

charles_npm.o :			$(CORE_DEPS) $(FLOW_HDR_DEPS) $(LSP_SOLVER_DEPS) \
				States.hpp BcCommon.hpp CtiMemory.hpp CommContainers.hpp \
				nonpremixed/NonpremixedSolver.hpp charles_npm.cpp

charles_npm.exe :		$(CHARLES_NPM_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

#-----------------------------------------
# helmholtz solver
#-----------------------------------------

HELMHOLTZ_OBJECTS=			helmholtz/HelmholtzSolver.o helmholtz/HelmholtzSolverSgs.o \
					helmholtz/HelmholtzSolverIo.o helmholtz/HelmholtzSolverBcs.o helmholtz/HelmholtzSolverWm.o

helmholtz/HelmholtzSolver.o:		helmholtz/HelmholtzSolver.hpp $(CORE_DEPS) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolver.cpp

helmholtz/HelmholtzSolverSgs.o:		helmholtz/HelmholtzSolver.hpp $(CORE_DEPS) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolverSgs.cpp

helmholtz/HelmholtzSolverWm.o:		helmholtz/HelmholtzSolver.hpp $(CORE_DEPS) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolverWm.cpp

helmholtz/HelmholtzSolverIo.o:		helmholtz/HelmholtzSolver.hpp $(CORE_DEPS) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolverIo.cpp

helmholtz/HelmholtzSolverBcs.o:		helmholtz/HelmholtzSolverBcs.hpp helmholtz/HelmholtzSolver.hpp bcprofile/ProfileReader.hpp \
					$(CORE_DEPS) $(FLOW_HDR_DEPS) helmholtz/HelmholtzSolverBcs.cpp

charles_h:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles_h..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles_h.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles_h.exe
	make ../../bin/charles_h.exe

CHARLES_H_OBJECTS = 		charles_h.o $(HELMHOLTZ_OBJECTS)

charles_h.o :			$(CORE_DEPS) $(FLOW_HDR_DEPS) $(LSP_SOLVER_DEPS) \
				helmholtz/HelmholtzSolver.hpp charles_h.cpp

charles_h.exe :		$(CHARLES_H_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

#---------------------------------
# main solver driver
#---------------------------------

CHARLES_OBJECTS = 		charles.o $(IG_OBJECTS) $(PREMIXED_OBJECTS) $(NONPREMIXED_OBJECTS) $(HELMHOLTZ_OBJECTS)

charles.o :			../core/common/libcti_core.a $(CORE_DEPS) $(FLOW_HDR_DEPS) \
				States.hpp BcCommon.hpp BasicPostpro.hpp CtiMemory.hpp CommContainers.hpp \
				ideal_gas/IdealGasSolver.hpp premixed/PremixedSolver.hpp nonpremixed/NonpremixedSolver.hpp \
				helmholtz/HelmholtzSolver.hpp $(LSP_SOLVER_DEPS) charles.cpp

charles.exe :			$(CHARLES_OBJECTS)
	$(CXX) -o $@ $^ $(CTI_LIBS)

# ==================================------===============================
# serial charles...

CXXFLAGS_NO_MPI	+= 	-I../core/common
CXXFLAGS_NO_MPI	+= 	-I../core/encrypt
CXXFLAGS_NO_MPI	+= 	-I../chemistry
CXXFLAGS_NO_MPI	+= 	-I../quiver
CXXFLAGS_NO_MPI	+=	-I./lsp
CXXFLAGS_NO_MPI	+= 	-I./ideal_gas
CXXFLAGS_NO_MPI	+= 	-I./helmholtz
CXXFLAGS_NO_MPI	+= 	-I./premixed
CXXFLAGS_NO_MPI	+= 	-I./nonpremixed
CXXFLAGS_NO_MPI	+= 	-I.

CXXFLAGS_NO_MPI	+=	$(CHARLES_OPT_FLAGS)

CTI_LIBS_S = 		-L../core/common -lcti_core_s -L./lsp -lcti_lsp_s \
			-L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

#Rule for serial compilation
%_s.o : %.cpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@
	@#@echo "     CXX_S $<"

charles_s:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for charles_s..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making charles_s.exe..."
	@echo "-------------------------------------------------"
	@echo
	make charles_s.exe

CORE_DEPS_S =			../core/common/libcti_core_s.a ../core/common/StaticSolver_write_data_s.o  \
				../core/common/StaticSolver_probes_s.o ../core/common/StaticSolver_multigrid_s.o \
				../core/common/StaticSolver_s.o ../core/common/Stl_s.o ../core/common/Multigrid_s.o

#LSP_SOLVER_DEPS_S=		lsp/IdealGasSolverWithLsp.hpp lsp/SingleLspSolver.hpp lsp/libcti_lsp_s.a
LSP_SOLVER_DEPS_S=		lsp/IdealGasSolverWithLsp.hpp lsp/libcti_lsp_s.a

lsp/libcti_lsp_s.a:       	$(LSP_HDR_DEPS)
	make -C ./lsp cti_lsp_s

#-----------------------------------------
# ideal gas solver
#-----------------------------------------

IG_OBJECTS_S=			ideal_gas/IdealGasSolver_s.o ideal_gas/IdealGasSolverSgs_s.o \
				ideal_gas/IdealGasSolverIo_s.o ideal_gas/IdealGasSolverBcs_s.o

ideal_gas/IdealGasSolver_s.o:		ideal_gas/IdealGasSolver.hpp ideal_gas/IdealGasSolverFlux.hpp \
					$(CORE_DEPS_S) $(FLOW_HDR_DEPS) ideal_gas/IdealGasSolver.cpp

ideal_gas/IdealGasSolverSgs_s.o:	ideal_gas/IdealGasSolver.hpp Sgs.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ideal_gas/IdealGasSolverSgs.cpp

ideal_gas/IdealGasSolverIo_s.o:		ideal_gas/IdealGasSolver.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ideal_gas/IdealGasSolverIo.cpp

ideal_gas/IdealGasSolverBcs_s.o:	ideal_gas/IdealGasSolverBcs.hpp ideal_gas/IdealGasSolver.hpp \
					ideal_gas/IdealGasSolverFlux.hpp \
					it/InflowTurbulenceIG.hpp BcCommon.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ideal_gas/IdealGasSolverBcs.cpp

#-----------------------------------------
# premixed solver
#-----------------------------------------

PREMIXED_OBJECTS_S=			premixed/PremixedSolver_s.o premixed/PremixedSolverSgs_s.o \
					premixed/PremixedSolverIo_s.o premixed/PremixedSolverBcs_s.o

premixed/PremixedSolver_s.o:		premixed/PremixedSolver.hpp premixed/PremixedSolverFlux.hpp \
					$(CORE_DEPS_S) $(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp  premixed/PremixedSolver.cpp

premixed/PremixedSolverSgs_s.o:		premixed/PremixedSolver.hpp Sgs.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp premixed/PremixedSolverSgs.cpp

premixed/PremixedSolverIo_s.o:		premixed/PremixedSolver.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp premixed/PremixedSolverIo.cpp

premixed/PremixedSolverBcs_s.o:		premixed/PremixedSolverBcs.hpp premixed/PremixedSolver.hpp \
					premixed/PremixedSolverFlux.hpp it/InflowTurbulenceP.hpp BcCommon.hpp \
					$(CORE_DEPS_S) $(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp premixed/PremixedSolverBcs.cpp

#-----------------------------------------
# nonpremixed solver
#-----------------------------------------

NONPREMIXED_OBJECTS_S=			nonpremixed/NonpremixedSolver_s.o nonpremixed/NonpremixedSolverSgs_s.o \
					nonpremixed/NonpremixedSolverIo_s.o nonpremixed/NonpremixedSolverBcs_s.o

nonpremixed/NonpremixedSolver_s.o:	nonpremixed/NonpremixedSolver.hpp premixed/PremixedSolverFlux.hpp \
					$(CORE_DEPS_S) $(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp  nonpremixed/NonpremixedSolver.cpp

nonpremixed/NonpremixedSolverSgs_s.o:	nonpremixed/NonpremixedSolver.hpp Sgs.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolverSgs.cpp

nonpremixed/NonpremixedSolverIo_s.o:	nonpremixed/NonpremixedSolver.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolverIo.cpp

nonpremixed/NonpremixedSolverBcs_s.o:	nonpremixed/NonpremixedSolverBcs.hpp nonpremixed/NonpremixedSolver.hpp \
					premixed/PremixedSolverFlux.hpp BcCommon.hpp $(CORE_DEPS_S) \
					$(FLOW_HDR_DEPS) ../chemistry/Chemtable.hpp nonpremixed/NonpremixedSolverBcs.cpp

#-----------------------------------------
# helmholtz solver
#-----------------------------------------

HELMHOLTZ_OBJECTS_S=			helmholtz/HelmholtzSolver_s.o helmholtz/HelmholtzSolverSgs_s.o \
					helmholtz/HelmholtzSolverIo_s.o helmholtz/HelmholtzSolverBcs_s.o helmholtz/HelmholtzSolverWm_s.o

helmholtz/HelmholtzSolver_s.o:		helmholtz/HelmholtzSolver.hpp $(CORE_DEPS_S) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolver.cpp

helmholtz/HelmholtzSolverSgs_s.o:	helmholtz/HelmholtzSolver.hpp $(CORE_DEPS_S) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolverSgs.cpp

helmholtz/HelmholtzSolverWm_s.o:	helmholtz/HelmholtzSolver.hpp $(CORE_DEPS_S) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolverWm.cpp

helmholtz/HelmholtzSolverIo_s.o:	helmholtz/HelmholtzSolver.hpp $(CORE_DEPS_S) $(FLOW_HDR_DEPS) \
					helmholtz/HelmholtzSolverIo.cpp

helmholtz/HelmholtzSolverBcs_s.o:	helmholtz/HelmholtzSolverBcs.hpp helmholtz/HelmholtzSolver.hpp \
					$(CORE_DEPS_S) $(FLOW_HDR_DEPS) helmholtz/HelmholtzSolverBcs.cpp

#---------------------------------
# main solver driver
#---------------------------------

CHARLES_OBJECTS_S = 		charles_s.o $(IG_OBJECTS_S) $(PREMIXED_OBJECTS_S) $(NONPREMIXED_OBJECTS_S) \
				$(HELMHOLTZ_OBJECTS_S) $(CORE_DEPS_S)

charles_s.o :		 	$(CORE_DEPS_S) $(FLOW_HDR_DEPS) States.hpp BcCommon.hpp BasicPostpro.hpp \
				CtiMemory.hpp CommContainers.hpp ideal_gas/IdealGasSolver.hpp premixed/PremixedSolver.hpp \
				nonpremixed/NonpremixedSolver.hpp helmholtz/HelmholtzSolver.hpp $(LSP_SOLVER_DEPS_S) charles.cpp

charles_s.exe :			$(CHARLES_OBJECTS_S)
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS_S)

# ==================================------===============================

testFlowSolver:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for testFlowSolver..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making testFlowSolver.exe..."
	@echo "-------------------------------------------------"
	@echo
	make testFlowSolver.exe

testFlowSolver.o :	$(CORE_DEPS) FlowSolver.hpp testFlowSolver.cpp

CTI_LIBS_NO_LSP =	-L../core/common -lcti_core \
			-L../core/encrypt -ltomlight $(CLIBS)

testFlowSolver.exe :	testFlowSolver.o
	$(CXX) -o $@ $^ $(CTI_LIBS_NO_LSP)

# ==================================------===============================

fea:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for fea..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core.a
	make -C ../core/encrypt
	@echo
	@echo "-------------------------------------------------"
	@echo "making fea.exe..."
	@echo "-------------------------------------------------"
	@echo
	make fea.exe
	make ../../bin/fea.exe

fea.o:	fea.cpp
fea.exe :	fea.o
	$(CXX) -o $@ $^ $(CTI_LIBS_NO_LSP)

# ==================================------===============================

tests:
	make -C tests

clean:
	make -C ./lsp clean
	rm -vf *.o
	rm -vf $(IG_OBJECTS) $(PREMIXED_OBJECTS) $(NONPREMIXED_OBJECTS) $(HELMHOLTZ_OBJECTS)
	rm -vf $(IG_OBJECTS_S) $(PREMIXED_OBJECTS_S) $(NONPREMIXED_OBJECTS) $(HELMHOLTZ_OBJECTS_S)
	rm -vf charles.exe ../../bin/charles.exe
	rm -vf charles_s.exe
	rm -vf charles_ig.exe
	rm -vf charles_pm.exe
	rm -vf charles_npm.exe
	rm -vf charles_h.exe
	rm -vf testFlowSolver.exe
	rm -vf tetMshToTetBin.exe
	rm -vf fea.exe

.PHONY: charles testFlowSolver tetMshToTetBin fea
