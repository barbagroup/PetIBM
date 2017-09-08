noinst_LTLIBRARIES = libutilities.la

libutilities_la_SOURCES = bodypack.cpp \
                          boundary.cpp \
                          cartesianmesh.cpp \
                          delta.cpp \
                          flowdescription.cpp \
                          parser.cpp \
                          simulationparameters.cpp \
                          singlebody.cpp \
                          singleboundary.cpp \
                          solutions.cpp \
                          timeintegration.cpp \
                          types.cpp

libutilities_la_CPPFLAGS = -I$(top_srcdir)/include \
                           -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                           -I$(top_srcdir)/external/boost-1.57.0

if WITH_AMGXWRAPPER
libutilities_la_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta2/src \
                            -I$(CUDA_DIR)/include \
                            -I$(AMGX_DIR)/examples \
                            -D HAVE_AMGX
endif
