noinst_LIBRARIES = libsolvers.a

libsolvers_a_SOURCES = createSolver.cpp \
                       navierStokes/NavierStokesSolver.cpp \
                       tairaColonius/TairaColoniusSolver.cpp

libsolvers_a_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                        -I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
                        -I$(top_srcdir)/external/boost-1.57.0