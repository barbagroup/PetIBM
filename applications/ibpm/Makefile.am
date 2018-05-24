bin_PROGRAMS = petibm-ibpm

petibm_ibpm_SOURCES = \
	main.cpp \
	ibpm.cpp

petibm_ibpm_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

petibm_ibpm_LDADD = \
	$(top_builddir)/applications/navierstokes/petibm_navierstokes-navierstokes.o \
	$(top_builddir)/src/libpetibm.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
