SUBDIRS = utilities \
          operators \
          linSolvers

lib_LTLIBRARIES = libpetibm.la

libpetibm_la_SOURCES =

libpetibm_la_LIBADD = $(top_builddir)/src/linSolvers/liblinsolvers.la \
                      $(top_builddir)/src/operators/liboperators.la \
                      $(top_builddir)/src/utilities/libutilities.la \
                      $(top_builddir)/external/yaml-cpp-0.5.1/libyamlcpp.la

if WITH_AMGXWRAPPER
libpetibm_la_LIBADD += $(top_builddir)/external/AmgXWrapper-1.0-beta2/libamgxwrapper.la
endif
