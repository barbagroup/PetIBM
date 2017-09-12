noinst_LTLIBRARIES = libyamlcpp.la

libyamlcpp_la_SOURCES = src/binary.cpp \
                        src/convert.cpp \
                        src/directives.cpp \
                        src/emit.cpp \
                        src/emitfromevents.cpp \
                        src/emitter.cpp \
                        src/emitterstate.cpp \
                        src/emitterutils.cpp \
                        src/exp.cpp \
                        src/memory.cpp \
                        src/node.cpp \
                        src/node_data.cpp \
                        src/nodebuilder.cpp \
                        src/nodeevents.cpp \
                        src/null.cpp \
                        src/ostream_wrapper.cpp \
                        src/parse.cpp \
                        src/parser.cpp \
                        src/regex.cpp \
                        src/scanner.cpp \
                        src/scanscalar.cpp \
                        src/scantag.cpp \
                        src/scantoken.cpp \
                        src/simplekey.cpp \
                        src/singledocparser.cpp \
                        src/stream.cpp \
                        src/tag.cpp \
                        src/contrib/graphbuilder.cpp \
                        src/contrib/graphbuilderadapter.cpp

libyamlcpp_la_CPPFLAGS = -I$(srcdir)/include \
                         -I$(srcdir)/src \
                         -I$(top_srcdir)/external/boost-1.57.0

install:
	mkdir -p $(includedir)/yaml-cpp
	cp -r $(top_srcdir)/external/yaml-cpp-0.5.1/include/yaml-cpp/* $(includedir)/yaml-cpp ;
	mkdir -p $(includedir)/boost
	cp -r $(top_srcdir)/external/boost-1.57.0/boost/* $(includedir)/boost ;
