if WITH_AMGXWRAPPER
noinst_LIBRARIES = libamgxwrapper.a
endif

libamgxwrapper_a_SOURCES = src/AmgXSolver.cpp

libamgxwrapper_a_CPPFLAGS = -I$(srcdir)/src \
                            -I$(AMGX_DIR)/examples \
                            -I$(CUDA_DIR)/include
