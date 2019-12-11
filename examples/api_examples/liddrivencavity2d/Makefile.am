EXTRA_PROGRAMS = cavity

cavity_SOURCES = \
	main.cpp

cavity_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS) \
	$(YAMLCPP_CPPFLAGS)

cavity_LDADD = \
	$(top_builddir)/src/libpetibm.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS) \
	$(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
if WITH_AMGX
cavity_LDADD += $(AMGXWRAPPER_LDFLAGS) $(AMGXWRAPPER_LIBS)
endif

cavity_DEPENDENCIES = input_data

input_data:
	if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		cp -f $(srcdir)/README.md $(PWD) ; \
		cp -f $(srcdir)/config.yaml $(PWD) ; \
		cp -rf $(srcdir)/config $(PWD) ; \
		cp -f $(srcdir)/postprocessing.xmf $(PWD) ; \
	fi ;

.PHONY: input_data
