bin_PROGRAMS = petibm-navierstokes

petibm_navierstokes_SOURCES = \
	main.cpp \
	navierstokes.cpp

petibm_navierstokes_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

petibm_navierstokes_LDADD = \
	$(top_builddir)/src/libpetibm.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
