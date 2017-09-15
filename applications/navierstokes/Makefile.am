bin_PROGRAMS = petibm-navierstokes

petibm_navierstokes_SOURCES = main.cpp \
                              navierstokes.cpp

petibm_navierstokes_CPPFLAGS = -I$(top_srcdir)/include \
                               $(YAMLCPP_CPPFLAGS)

petibm_navierstokes_LDADD = $(top_builddir)/src/libpetibm.la \
                            $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
petibm_navierstokes_CPPFLAGS += $(AMGXWRAPPER_CPPFLAGS)
petibm_navierstokes_LDADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif

