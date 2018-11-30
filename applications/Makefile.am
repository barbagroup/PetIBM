SUBDIRS = \
	navierstokes \
	ibpm \
	decoupledibpm \
	vorticity \
	createxdmf \
	writemesh

lib_LTLIBRARIES = libpetibmapps.la

libpetibmapps_la_SOURCES = \
	navierstokes/navierstokes.cpp \
	ibpm/ibpm.cpp \
	decoupledibpm/decoupledibpm.cpp \
	rigidkinematics/rigidkinematics.cpp

libpetibmapps_la_CPPFLAGS = \
	-I$(top_srcdir)/include \
	navierstokes \
	ibpm \
	decoupledibpm \
	rigidkinematics \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

libpetibmapps_la_LIBADD = \
	$(top_builddir)/src/libpetibm.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)

if WITH_AMGX
libpetibmapps_la_LIBADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif

navierstokesincludedir = $(includedir)/petibm/navierstokes
navierstokesinclude_HEADERS = navierstokes/navierstokes.h

ibpmincludedir = $(includedir)/petibm/ibpm
ibpminclude_HEADERS = ibpm/ibpm.h

decoupledibpmincludedir = $(includedir)/petibm/decoupledibpm
decoupledibpminclude_HEADERS = decoupledibpm/decoupledibpm.h

rigidkinematicsincludedir = $(includedir)/petibm/rigidkinematics
rigidkinematicsinclude_HEADERS = rigidkinematics/rigidkinematics.h
