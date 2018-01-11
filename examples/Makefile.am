copy-examples:
	@echo "Copying examples folders..."
	@if test "$(top_srcdir)" != "$(top_builddir)" ; then \
		rsync -a $(top_srcdir)/examples/* $(PWD) \
			--exclude decoupledibpm/oscillatingcylinder2dRe100_GPU ; \
	fi ;

.PHONY: input_data
