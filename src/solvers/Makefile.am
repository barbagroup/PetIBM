noinst_LIBRARIES = libsolvers.a
libsolvers_a_SOURCES = createSolver.cpp \
                       navierStokes/NavierStokesSolver.cpp \
                       tairaColonius/TairaColoniusSolver.cpp
libsolvers_a_CPPFLAGS = -I$(top_srcdir)/src/utilities \
                        -I$(top_srcdir)/external/yaml-cpp-0.5.1/include
libsolvers_a_CXXFLAGS = -std=c++0x