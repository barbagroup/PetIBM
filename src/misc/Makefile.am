noinst_LTLIBRARIES = libmisc.la

libmisc_la_SOURCES = \
	misc.cpp \
	type.cpp \
	delta.cpp

libmisc_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libmisc_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
