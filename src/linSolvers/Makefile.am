noinst_LIBRARIES = liblinsolvers.a

liblinsolvers_a_SOURCES = KSPSolver.cpp \
												  LinSolver.cpp

liblinsolvers_a_CPPFLAGS = -I$(top_srcdir)/src/utilities

if WITH_AMGXWRAPPER
liblinsolvers_a_SOURCES += AMGXSolver.cpp
liblinsolvers_a_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta2/src \
                            -I$(CUDA_DIR)/include \
                            -I$(AMGX_DIR)/examples \
                            -D HAVE_AMGX
endif
