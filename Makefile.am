ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTION = foreign

SUBDIRS = \
	src \
	tests \
	applications \
	include

EXTRA_DIST = autogen.sh

EXAMPLES_DIR=$(top_builddir)/examples

copy-examples:
	@if test $(EXAMPLES_DIR) != "$(top_srcdir)/examples" ; then \
		echo "Copying examples folders into $(EXAMPLES_DIR)..." ; \
		mkdir -p $(EXAMPLES_DIR) ; \
		rsync -a $(top_srcdir)/examples/* $(EXAMPLES_DIR) ; \
	fi ;

.PHONY: copy-examples
