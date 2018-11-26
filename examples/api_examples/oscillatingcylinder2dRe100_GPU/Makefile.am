EXTRA_PROGRAMS = oscillatingcylinder

oscillatingcylinder_SOURCES = \
	main.cpp \
	oscillatingcylinder.cpp

oscillatingcylinder_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

oscillatingcylinder_LDADD = \
	$(top_builddir)/applications/libpetibmapps.la \
	$(top_builddir)/src/libpetibm.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
if WITH_AMGX
oscillatingcylinder_LDADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif

oscillatingcylinder_DEPENDENCIES = input_data

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -f $(srcdir)/README.md $(PWD) ; \
		cp -f $(srcdir)/config.yaml $(PWD) ; \
		cp -rf $(srcdir)/config $(PWD) ; \
		cp -f $(srcdir)/circle.body $(PWD) ; \
		cp -rf $(srcdir)/scripts $(PWD) ; \
	fi ;

.PHONY: input_data
