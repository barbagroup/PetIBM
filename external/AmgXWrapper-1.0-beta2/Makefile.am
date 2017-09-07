if WITH_AMGXWRAPPER
noinst_LTLIBRARIES = libamgxwrapper.la
endif

libamgxwrapper_la_SOURCES = src/AmgXSolver.cpp

libamgxwrapper_la_CPPFLAGS = -I$(srcdir)/src \
                             -I$(AMGX_DIR)/examples \
                             -I$(CUDA_DIR)/include
