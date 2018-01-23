bin_PROGRAMS = petibm-vorticity

petibm_vorticity_SOURCES = main.cpp

petibm_vorticity_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

petibm_vorticity_LDADD = \
	$(top_builddir)/src/libpetibm.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
