noinst_LIBRARIES = libutilities.a

libutilities_a_SOURCES = CartesianMesh.cpp \
                         FlowDescription.cpp \
                         SimulationParameters.cpp \
                         Body.cpp \
                         types.cpp \
                         solvers/kspsolver.cpp

if WITH_AMGXWRAPPER
libutilities_a_SOURCES += solvers/amgxsolver.cpp
endif

libutilities_a_CPPFLAGS = -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                          -I$(top_srcdir)/external/boost-1.57.0


if WITH_AMGXWRAPPER
libutilities_a_CPPFLAGS += -I$(top_srcdir)/external/AmgXWrapper-1.0-beta/src \
                           -I$(CUDA_DIR)/include \
                           -I$(AMGX_DIR)/examples
endif
