bin_PROGRAMS = petibm-createxdmf

petibm_createxdmf_SOURCES = main.cpp

petibm_createxdmf_CPPFLAGS = -I$(top_srcdir)/include \
                             $(PETSC_CPPFLAGS) \
                             $(YAMLCPP_CPPFLAGS)

petibm_createxdmf_LDADD = $(top_builddir)/src/libpetibm.la \
                         $(PETSC_LDFLAGS) $(PETSC_LIBS) \
                         $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
