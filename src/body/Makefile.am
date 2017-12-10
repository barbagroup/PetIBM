noinst_LTLIBRARIES = libbody.la

libbody_la_SOURCES = bodypack.cpp singlebody.cpp singlebodypoints.cpp

libbody_la_CPPFLAGS = -I$(top_srcdir)/include $(PETSC_CPPFLAGS) $(YAMLCPP_CPPFLAGS)

libbody_la_LIBADD = $(PETSC_LDFLAGS) $(PETSC_LIBS) $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
