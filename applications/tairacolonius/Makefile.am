bin_PROGRAMS = petibm-tairacolonius

petibm_tairacolonius_SOURCES = main.cpp \
                               tairacolonius.cpp

petibm_tairacolonius_CPPFLAGS = -I$(top_srcdir)/include \
                               $(PETSC_CPPFLAGS) \
                               $(YAMLCPP_CPPFLAGS)

petibm_tairacolonius_LDADD = $(top_builddir)/applications/navierstokes/petibm_navierstokes-navierstokes.o \
                            $(top_builddir)/src/libpetibm.la \
                            $(PETSC_LDFLAGS) $(PETSC_LIBS) \
                            $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
