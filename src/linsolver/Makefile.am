noinst_LTLIBRARIES = liblinsolver.la

liblinsolver_la_SOURCES = \
	linsolver.cpp \
	linsolverksp.cpp

liblinsolver_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) $(YAMLCPP_CPPFLAGS)

liblinsolver_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
liblinsolver_la_SOURCES += linsolveramgx.cpp
liblinsolver_la_CPPFLAGS += -DHAVE_AMGX $(AMGXWRAPPER_CPPFLAGS)
liblinsolver_la_LIBADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif
