bin_PROGRAMS = petibm-decoupledibpm

petibm_decoupledibpm_SOURCES = main.cpp \
                               decoupledibpm.cpp

petibm_decoupledibpm_CPPFLAGS = -I$(top_srcdir)/include \
                                $(YAMLCPP_CPPFLAGS)

petibm_decoupledibpm_LDADD = $(top_builddir)/src/libpetibm.la \
                             $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
petibm_decoupledibpm_CPPFLAGS += $(AMGXWRAPPER_CPPFLAGS)
petibm_decoupledibpm_LDADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif

