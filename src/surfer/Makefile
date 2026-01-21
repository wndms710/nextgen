include ../../Makefile.in

# includes for the core (needs to be updated so its the include dir)
CXXFLAGS_NO_MPI +=	-I../core/common
CXXFLAGS_NO_MPI +=	-I../core/encrypt

CXXFLAGS +=	-I../core/common
CXXFLAGS +=	-I../core/encrypt

CTI_LIBS = 	-L../core/common -lcti_core_s -L../core/encrypt -ltomlight_s $(CLIBS_NO_MPI)

default: 	surfer post


../../bin/%.exe : %.exe
	cp $^ ../../bin/

surfer:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for surfer..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making surfer.exe..."
	@echo "-------------------------------------------------"
	@echo
	make surfer.exe
	make ../../bin/surfer.exe

post:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for post..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making post.exe..."
	@echo "-------------------------------------------------"
	@echo
	make post.exe
	make ../../bin/post.exe


test:
	@echo
	@echo "-------------------------------------------------"
	@echo "making components for test_surfer..."
	@echo "-------------------------------------------------"
	@echo
	make -C ../core/common libcti_core_s.a
	make -C ../core/encrypt libtomlight_s.a
	@echo
	@echo "-------------------------------------------------"
	@echo "making test_surfer.exe..."
	@echo "-------------------------------------------------"
	@echo
	make test_surfer.exe



diffgeom:
	@echo
	@echo "-------------------------------------------------"
	@echo "making diffgeom.exe..."
	@echo "-------------------------------------------------"
	@echo
	make diffgeom.exe

#Rule for serial compilation
%_s.o : %.cpp
	$(CXX_NO_MPI) -DNO_MPI $(CXXFLAGS_NO_MPI) $< -c -o $@
	@#@echo "     CXX_S $<"


../core/common/libcti_core_s.a:
	make -C ../core/common libcti_core_s.a

# ==================================
# SimpleSurface class
# ==================================

DirectedGraph_s.o:	../core/common/libcti_core_s.a DirectedGraph.hpp DirectedGraph.cpp

SimpleSurface_sbin_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_sbin.cpp

SimpleSurface_stl_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_stl.cpp DoubleVertex.hpp AsciiReader.hpp

SimpleSurface_ply_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_ply.cpp PlyReader.hpp

SimpleSurface_ugrid_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_ugrid.cpp UgridReader.hpp AsciiReader.hpp

SimpleSurface_obj_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_obj.cpp DoubleVertex.hpp

SimpleSurface_msh_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_msh.cpp ../core/common/parsers/fluent/CasReader.hpp ../core/common/parsers/fluent/FluentParser.hpp

SimpleSurface_restart_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_restart.cpp

SimpleSurface_primitives_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_primitives.cpp

SimpleSurface_lifted_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_lifted.cpp

SimpleSurface_naca_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_naca.cpp

SimpleSurface_tecplot_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_tecplot.cpp

SimpleSurface_input_files_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_input_files.cpp

SimpleSurface_data_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_data.cpp

SimpleSurface_periodic_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_periodic.cpp

SimpleSurface_transform_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_transform.cpp

SimpleSurface_geom_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_geom.cpp

SimpleSurface_edges_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_edges.cpp

SimpleSurface_teost_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_teost.cpp Teost.hpp DoubleVertex.hpp

SimpleSurface_imprint_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_imprint.cpp Teost.hpp DoubleVertex.hpp

SimpleSurface_repair_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp DirectedGraph.hpp SimpleSurface_repair.cpp

SimpleSurface_zipper_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_zipper.cpp

SimpleSurface_zoning_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_zoning.cpp

SimpleSurface_material_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_material.cpp

SimpleSurface_morph_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_morph.cpp

SimpleSurface_seed_s.o:	../core/common/libcti_core_s.a ../stitch/CuttableVoronoiData.hpp ../stitch/startEdge.hpp \
	../stitch/completeEdge.hpp ../stitch/checkOrSetEdge.hpp SimpleSurface.hpp SimpleSurface_seed.cpp

SimpleSurface_intersect_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_intersect.cpp

SimpleSurface_vor2d_s.o:	../core/common/libcti_core_s.a SimpleSurface.hpp SimpleSurface_vor2d.cpp

SIMPLE_SURFACE_OBJS_S = SimpleSurface_sbin_s.o SimpleSurface_stl_s.o SimpleSurface_msh_s.o SimpleSurface_ply_s.o SimpleSurface_ugrid_s.o SimpleSurface_obj_s.o \
	SimpleSurface_restart_s.o SimpleSurface_primitives_s.o SimpleSurface_lifted_s.o SimpleSurface_naca_s.o SimpleSurface_tecplot_s.o \
	SimpleSurface_input_files_s.o SimpleSurface_geom_s.o SimpleSurface_edges_s.o SimpleSurface_teost_s.o SimpleSurface_imprint_s.o \
	SimpleSurface_data_s.o SimpleSurface_periodic_s.o SimpleSurface_transform_s.o SimpleSurface_repair_s.o SimpleSurface_zipper_s.o \
	SimpleSurface_zoning_s.o SimpleSurface_material_s.o SimpleSurface_morph_s.o SimpleSurface_seed_s.o DirectedGraph_s.o \
	SimpleSurface_intersect_s.o SimpleSurface_vor2d_s.o

# ==================================
# SimpleRestart class
# ==================================

SimpleRestart_s.o:	../core/common/libcti_core_s.a SimpleRestart.hpp SimpleRestart.cpp


# ==================================
# surfer
# ==================================

surfer_s.o :	../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) NonManifoldTriFlag.hpp  testSeed.hpp surfer.siData surfer.cpp

surfer.exe :	../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) surfer_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)


# ==================================
# post
# ==================================

post_s.o :	../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) SimpleRestart_s.o post.siData post.cpp

post.exe :	../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) SimpleRestart_s.o post_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)



# ==================================
# diffgeom
# ==================================

diffgeom_s.o:	../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) diffgeom.cpp

diffgeom.exe : 	../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) diffgeom_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)

# ==================================
# testMinHeap
# ==================================

testMinHeap_s.o:		../core/common/libcti_core_s.a MinHeap.hpp testMinHeap.cpp

testMinHeap.exe:	../core/common/libcti_core_s.a testMinHeap_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)

# ==================================
# testDirectedGraph
# ==================================

testDirectedGraph_s.o:		../core/common/libcti_core_s.a DirectedGraph.hpp testDirectedGraph.cpp

testDirectedGraph.exe:	../core/common/libcti_core_s.a testDirectedGraph_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)


# ==================================
# test_surfer
# ==================================

tests/testDirectedGraph_s.o: tests/testDirectedGraph.cpp

tests/testSimpleSurface_s.o: $(SIMPLE_SURFACE_OBJS_S) tests/testSimpleSurface.cpp

TEST_OBJS= tests/testSimpleSurface_s.o \
	   tests/testDirectedGraph_s.o

test/test_surfer_s.o : ../core/catch1/catch.hpp tests/test_surfer.cpp

test_surfer.exe : ../core/common/libcti_core_s.a $(SIMPLE_SURFACE_OBJS_S) $(TEST_OBJS)  tests/test_surfer_s.o
	$(CXX_NO_MPI) -o $@ $^ $(CTI_LIBS)



clean:
	rm -vf *.o tests/*.o
	rm -vf surfer.exe test_surfer.exe diffgeom.exe testMinHeap.exe post.exe

vclean:
	make -C ../core/common clean
	make clean

.PHONY:	surfer post test diffgeom testMinHeap
