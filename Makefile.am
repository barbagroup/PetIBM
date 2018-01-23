ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTION = foreign

SUBDIRS = \
	src \
	tests \
	applications \
	include

EXTRA_DIST = autogen.sh

copy-examples:
	@echo "Copying examples folders..."
	@if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		mkdir -p $(PWD)/examples ;\
		rsync -a $(top_srcdir)/examples/* examples \
			--exclude decoupledibpm/oscillatingcylinder2dRe100_GPU ; \
	fi ;

.PHONY: copy-examples
