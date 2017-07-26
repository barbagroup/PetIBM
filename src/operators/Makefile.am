noinst_LIBRARIES = liboperators.a

liboperators_a_SOURCES = createBn.cpp \
												 createConvection.cpp \
												 createDelta.cpp \
												 createDiagMatrix.cpp \
                         createDivergence.cpp \
                         createGradient.cpp \
                         createLaplacian.cpp

liboperators_a_CPPFLAGS = -I$(top_srcdir)/src/utilities \
													-I$(top_srcdir)/external/yaml-cpp-0.5.1/include \
													-I$(top_srcdir)/external/boost-1.57.0
