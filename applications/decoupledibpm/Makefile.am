bin_PROGRAMS = petibm-decoupledibpm

petibm_decoupledibpm_SOURCES = main.cpp \
                               decoupledibpm.cpp

petibm_decoupledibpm_CPPFLAGS = -I$(top_srcdir)/include \
                                $(PETSC_CPPFLAGS) \
                                $(YAMLCPP_CPPFLAGS)

petibm_decoupledibpm_LDADD = $(top_builddir)/applications/navierstokes/petibm_navierstokes-navierstokes.o \
                             $(top_builddir)/src/libpetibm.la \
                             $(PETSC_LDFLAGS) $(PETSC_LIBS) \
                             $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
