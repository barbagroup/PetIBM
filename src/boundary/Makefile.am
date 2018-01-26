noinst_LTLIBRARIES = libboundary.la

libboundary_la_SOURCES = \
	boundary.cpp \
	boundarysimple.cpp \
	singleboundary.cpp \
	singleboundarybase.cpp \
	singleboundaryconvective.cpp \
	singleboundarydirichlet.cpp \
	singleboundaryneumann.cpp \
	singleboundaryperiodic.cpp

libboundary_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libboundary_la_LIBADD = \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
