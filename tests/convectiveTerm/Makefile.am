include $(PETSC_DIR)/conf/variables

bin_PROGRAMS = convectiveTermTest2d \
							 convectiveTermTest3d

convectiveTermTest2d_SOURCES = ConvectiveTerm.cpp \
															 convectiveTermTest.cpp

convectiveTermTest2d_DEPENDENCIES = input_data

convectiveTermTest2d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
																-I$(top_srcdir)/src/solvers \
																-D DIMENSIONS=2

convectiveTermTest2d_CXXFLAGS = -std=c++0x

convectiveTermTest2d_LDADD = $(top_builddir)/src/utilities/libutilities.a \
														 $(top_builddir)/src/solvers/libsolvers.a \
														 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

convectiveTermTest3d_SOURCES = ConvectiveTerm.cpp \
															 convectiveTermTest.cpp

convectiveTermTest3d_DEPENDENCIES = input_data

convectiveTermTest3d_CPPFLAGS = -I$(top_srcdir)/src/utilities \
																-I$(top_srcdir)/src/solvers \
																-D DIMENSIONS=3

convectiveTermTest3d_CXXFLAGS = -std=c++0x

convectiveTermTest3d_LDADD = $(top_builddir)/src/utilities/libutilities.a \
														 $(top_builddir)/src/solvers/libsolvers.a \
														 $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.a

input_data:
	@echo "Copying cases folder..."
	@if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -r $(top_srcdir)/tests/convectiveTerm/cases $(PWD) ; \
		cp $(top_srcdir)/tests/convectiveTerm/plotGridConvergence.py $(PWD) ; \
	fi ;

tests:
	@echo "Running 2D convective-term tests..."
	@mkdir -p  data
	@$(RM) data/relativeErrors2d.dat
	@$(MPIEXEC) -n 1 ./convectiveTermTest2d -caseFolder ./cases/2d/6 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@$(MPIEXEC) -n 1 ./convectiveTermTest2d -caseFolder ./cases/2d/12 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@$(MPIEXEC) -n 1 ./convectiveTermTest2d -caseFolder ./cases/2d/24 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@$(MPIEXEC) -n 1 ./convectiveTermTest2d -caseFolder ./cases/2d/48 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@echo "Running 3D convective-term tests..."
	@$(RM) -f data/relativeErrors3d.dat
	@$(MPIEXEC) -n 1 ./convectiveTermTest3d -caseFolder ./cases/3d/6 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@$(MPIEXEC) -n 1 ./convectiveTermTest3d -caseFolder ./cases/3d/12 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@$(MPIEXEC) -n 1 ./convectiveTermTest3d -caseFolder ./cases/3d/24 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@$(MPIEXEC) -n 1 ./convectiveTermTest3d -caseFolder ./cases/3d/48 \
							-sys2_pc_type gamg -sys2_pc_gamg_type agg -sys2_pc_gamg_agg_nsmooths 1
	@echo "Plotting grid convergence for convective-term tests..."
	@python plotGridConvergence.py

.PHONY: input_data tests