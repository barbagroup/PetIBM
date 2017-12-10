noinst_LTLIBRARIES = libmesh.la

libmesh_la_SOURCES = \
	mesh.cpp \
	cartesianmesh.cpp

libmesh_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libmesh_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
