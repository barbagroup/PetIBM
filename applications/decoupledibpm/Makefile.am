bin_PROGRAMS = petibm-decoupledibpm

petibm_decoupledibpm_SOURCES = main.cpp \
                               decoupledibpm.cpp

petibm_decoupledibpm_CPPFLAGS = -I$(top_srcdir)/src \
                                -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                                -I$(top_srcdir)/external/boost-1.57.0

petibm_decoupledibpm_LDADD = $(top_builddir)/src/libpetibm.la \
                             -lstdc++fs
