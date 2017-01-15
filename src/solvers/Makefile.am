noinst_LIBRARIES = libsolvers.a

libsolvers_a_SOURCES = createSolver.cpp \
                       navierStokes/NavierStokesSolver.cpp \
                       tairaColonius/TairaColoniusSolver.cpp

libsolvers_a_CPPFLAGS = -I$(top_srcdir)/src/utilities

if WITH_AMGXWRAPPER
libsolvers_a_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta/src \
                         -I$(CUDA_DIR)/include \
                         -I$(AMGX_DIR)/examples \
                         -D HAVE_AMGX
endif
