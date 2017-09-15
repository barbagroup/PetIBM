noinst_LTLIBRARIES = liblinsolvers.la

liblinsolvers_la_SOURCES = kspsolver.cpp \
                           linsolver.cpp

liblinsolvers_la_CPPFLAGS = -I$(top_srcdir)/include \
                            $(YAMLCPP_CPPFLAGS)

liblinsolvers_la_LIBADD = $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
liblinsolvers_la_SOURCES += amgxsolver.cpp
liblinsolvers_la_CPPFLAGS += $(AMGXWRAPPER_CPPFLAGS)
liblinsolvers_la_LIBADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif
