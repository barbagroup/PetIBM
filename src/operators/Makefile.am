noinst_LTLIBRARIES = liboperators.la

liboperators_la_SOURCES = createbn.cpp \
                          createconvection.cpp \
                          createdelta.cpp \
                          creatediagmatrix.cpp \
                          createdivergence.cpp \
                          creategradient.cpp \
                          createlaplacian.cpp

liboperators_la_CPPFLAGS = -I$(top_srcdir)/include \
                           $(YAMLCPP_CPPFLAGS)

liboperators_la_LIBADD = -lstdc++fs \
                         $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
