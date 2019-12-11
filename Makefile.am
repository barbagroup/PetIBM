ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTION = foreign

SUBDIRS = \
	src \
	tests \
	applications \
	include

EXTRA_DIST = autogen.sh

install-data-hook:
	@if test "$(prefix)/examples" != "$(top_srcdir)/examples" ; then \
		echo "Copying examples into $(prefix)..." ; \
		rsync -a $(top_srcdir)/examples $(prefix) ; \
		rm -rf $(prefix)/examples/api_examples ; \
	fi ;
