noinst_LTLIBRARIES = liblinsolvers.la

liblinsolvers_la_SOURCES = kspsolver.cpp \
                           linsolver.cpp

liblinsolvers_la_CPPFLAGS = -I$(top_srcdir)/include \
                            $(PETSC_CPPFLAGS)

liblinsolvers_la_LIBADD = $(PETSC_LDFLAGS) $(PETSC_LIBS)

if WITH_AMGX
liblinsolvers_la_SOURCES += amgxsolver.cpp
liblinsolvers_la_CPPFLAGS += -DHAVE_AMGX $(AMGXWRAPPER_CPPFLAGS)
liblinsolvers_la_LIBADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif
