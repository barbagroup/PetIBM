noinst_LTLIBRARIES = libtimeintegration.la

libtimeintegration_la_SOURCES = \
	timeintegration.cpp

libtimeintegration_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libtimeintegration_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
