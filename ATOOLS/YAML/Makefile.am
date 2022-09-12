SUBDIRS = yaml-cpp

pkglib_LTLIBRARIES = libToolsYaml.la

EXTRA_DIST = yaml-cpp

libToolsYaml_la_SOURCES = \
			  binary.cpp \
			  convert.cpp \
			  depthguard.cpp \
			  directives.cpp \
			  emit.cpp \
			  emitfromevents.cpp \
			  emitter.cpp \
			  emitterstate.cpp \
			  emitterutils.cpp \
			  exceptions.cpp \
			  exp.cpp \
			  memory.cpp \
			  node.cpp \
			  node_data.cpp \
			  nodebuilder.cpp \
			  nodeevents.cpp \
			  null.cpp \
			  ostream_wrapper.cpp \
			  parse.cpp \
			  parser.cpp \
			  regex_yaml.cpp \
			  scanner.cpp \
			  scanscalar.cpp \
			  scantag.cpp \
			  scantoken.cpp \
			  simplekey.cpp \
			  singledocparser.cpp \
			  stream.cpp \
			  tag.cpp

localinc_HEADERS = \
		   collectionstack.h \
		   directives.h \
		   emitterstate.h \
		   emitterutils.h \
		   exp.h \
		   indentation.h \
		   nodebuilder.h \
		   nodeevents.h \
		   ptr_vector.h \
		   regex_yaml.h \
		   regeximpl.h \
		   scanner.h \
		   scanscalar.h \
		   scantag.h \
		   setting.h \
		   singledocparser.h \
		   stream.h \
		   streamcharsource.h \
		   stringsource.h \
		   tag.h \
		   token.h

libToolsYaml_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(srcdir)
