SUBDIRS = utilities \
          operators \
          linsolvers

lib_LTLIBRARIES = libpetibm.la

libpetibm_la_SOURCES =

libpetibm_la_LIBADD = $(top_builddir)/src/linsolvers/liblinsolvers.la \
                      $(top_builddir)/src/operators/liboperators.la \
                      $(top_builddir)/src/utilities/libutilities.la \
                      $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
libpetibm_la_LIBADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif
