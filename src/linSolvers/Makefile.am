noinst_LTLIBRARIES = liblinsolvers.la

liblinsolvers_la_SOURCES = KSPSolver.cpp \
                           LinSolver.cpp

liblinsolvers_la_CPPFLAGS = -I$(top_srcdir)/src

if WITH_AMGXWRAPPER
liblinsolvers_la_SOURCES += AMGXSolver.cpp
liblinsolvers_la_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta2/src \
                             -I$(CUDA_DIR)/include \
                             -I$(AMGX_DIR)/examples \
                             -D HAVE_AMGX
endif
