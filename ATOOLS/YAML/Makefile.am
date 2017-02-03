## How to make our mangled yaml-cpp from the yaml-cpp-0.3.0.tar.gz tarball
 #tar xf yaml-cpp-0.3.0.tar.gz
 #for i in `find . -type f -name "*.h" -o -name "*.cpp"`; do sed -e 's/namespace YAML/namespace SHERPA_YAML/' -e 's/YAML::/SHERPA_YAML::/g' -i; done
 #for i in `find . -type f -name "*.h" -o -name "*.cpp"`; do sed -i -e 's/namespace YAML/namespace SHERPA_YAML/' -e 's/YAML::/SHERPA_YAML::/g' $i; done
 #mv yaml-cpp{,.bak} && cp yaml-cpp.bak/include/yaml-cpp . -r && cp yaml-cpp.bak/src/*.{h,cpp} .

pkglib_LTLIBRARIES = libToolsYaml.la

#EXTRA_DIST = yaml-cpp
#noinst_LTLIBRARIES = libsherpa-yaml-cpp.la
libToolsYaml_la_SOURCES = \
  binary.cpp      emit.cpp            emitterstate.cpp  memory.cpp       \
  convert.cpp     emitfromevents.cpp  emitterutils.cpp  nodebuilder.cpp  \
  directives.cpp  emitter.cpp         exp.cpp           node.cpp        \
  node_data.cpp   ostream_wrapper.cpp  regex_yaml.cpp  scantag.cpp    singledocparser.cpp \
  nodeevents.cpp  parse.cpp            scanner.cpp     scantoken.cpp  stream.cpp \
  null.cpp        parser.cpp           scanscalar.cpp  simplekey.cpp  tag.cpp

localinc_HEADERS = \
  collectionstack.h  emitterstate.h  exp.h          nodebuilder.h \
  ptr_stack.h   regeximpl.h   scanner.h     scantag.h  singledocparser.h \
  stream.h        tag.h directives.h       emitterutils.h  indentation.h  \
  nodeevents.h   ptr_vector.h  regex_yaml.h  scanscalar.h  setting.h  \
  streamcharsource.h  stringsource.h  token.h

libToolsYaml_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir) -Wno-deprecated-declarations
