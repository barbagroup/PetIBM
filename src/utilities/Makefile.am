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
                           $(PETSC_CPPFLAGS) \
                           $(YAMLCPP_CPPFLAGS)
libutilities_la_LIBADD = $(PETSC_LDFLAGS) $(PETSC_LIBS) \
                         $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS) \
                         -lstdc++fs

if WITH_AMGX
libutilities_la_CPPFLAGS += -DHAVE_AMGX
endif
