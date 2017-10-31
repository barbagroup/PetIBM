noinst_LTLIBRARIES = libsolution.la

libsolution_la_SOURCES = \
	solution.cpp \
	solutionsimple.cpp

libsolution_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libsolution_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
