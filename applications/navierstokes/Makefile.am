bin_PROGRAMS = petibm-navierstokes

petibm_navierstokes_SOURCES = main.cpp \
                              navierstokes.cpp

petibm_navierstokes_CPPFLAGS = -I$(top_srcdir)/include \
                               -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                               -I$(top_srcdir)/external/boost-1.57.0

petibm_navierstokes_LDADD = $(top_builddir)/src/libpetibm.la \
                            -lstdc++fs
