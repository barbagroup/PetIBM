SUBDIRS = \
	body \
	boundary \
	io \
	linsolver \
	mesh \
	misc \
	operators \
	parser \
	solution \
	timeintegration

lib_LTLIBRARIES = libpetibm.la

libpetibm_la_SOURCES =

libpetibm_la_LIBADD = \
	$(top_builddir)/src/body/libbody.la \
	$(top_builddir)/src/boundary/libboundary.la \
	$(top_builddir)/src/io/libio.la \
	$(top_builddir)/src/linsolver/liblinsolver.la \
	$(top_builddir)/src/mesh/libmesh.la \
	$(top_builddir)/src/misc/libmisc.la \
	$(top_builddir)/src/operators/liboperators.la \
	$(top_builddir)/src/parser/libparser.la \
	$(top_builddir)/src/solution/libsolution.la \
	$(top_builddir)/src/timeintegration/libtimeintegration.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
libpetibm_la_LIBADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif
