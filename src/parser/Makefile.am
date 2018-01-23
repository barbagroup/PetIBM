noinst_LTLIBRARIES = libparser.la

libparser_la_SOURCES = \
	parser.cpp \
	yaml.cpp

libparser_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libparser_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
