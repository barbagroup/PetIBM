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
                           $(YAMLCPP_CPPFLAGS)
libutilities_la_LIBADD = -lstdc++fs \
                         $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
