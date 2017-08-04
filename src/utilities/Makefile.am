noinst_LTLIBRARIES = libutilities.la

libutilities_la_SOURCES = BodyPack.cpp \
                          Boundary.cpp \
                          CartesianMesh.cpp \
                          delta.cpp \
                          FlowDescription.cpp \
                          parser.cpp \
                          SimulationParameters.cpp \
                          SingleBody.cpp \
                          SingleBoundary.cpp \
                          Solutions.cpp \
                          TimeIntegration.cpp \
                          types.cpp

libutilities_la_CPPFLAGS = -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                           -I$(top_srcdir)/external/boost-1.57.0

if WITH_AMGXWRAPPER
libutilities_la_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta2/src \
                            -I$(CUDA_DIR)/include \
                            -I$(AMGX_DIR)/examples \
                            -D HAVE_AMGX
endif
