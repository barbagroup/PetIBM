noinst_LTLIBRARIES = liboperators.la

liboperators_la_SOURCES = createbn.cpp \
                          createconvection.cpp \
                          createdelta.cpp \
                          creatediagmatrix.cpp \
                          createdivergence.cpp \
                          creategradient.cpp \
                          createlaplacian.cpp

liboperators_la_CPPFLAGS = -I$(top_srcdir)/include \
                           $(PETSC_CPPFLAGS) \
                           $(YAMLCPP_CPPFLAGS)

liboperators_la_LIBADD = $(PETSC_LDFLAGS) $(PETSC_LIBS) \
                         $(YAMLCPP_LDFLAGS) $(YAMLCPP_LIBS)
