
SVNTAG = AddOns/Python
include ../../svn.make

noinst_HEADERS = \
	Sherpa.i \
	Initialization_Handler.i \
	Matrix_Element_Handler.i \
	Process_Base.i \
	Cluster_Amplitude.i \
	Flavour.i \
	Vec4.i \
	Exception.i \
	Message.i \
	iostream.i

if ENABLE_PYEXT

Sherpa_wrap.cxx: $(noinst_HEADERS)
	@$(SWIG) -c++ -python -I$(top_srcdir) $<

pkgpython_LTLIBRARIES = _Sherpa.la

_Sherpa_la_SOURCES = SVN_Info.C Sherpa_wrap.cxx
_Sherpa_la_CPPFLAGS= @PYTHON_CPPFLAGS@
_Sherpa_la_LDFLAGS = -module
_Sherpa_la_LIBADD  = -L@SHERPABUILDDIR@/Main -lSherpaMain \
		$(HZTOOLLIBS) $(PGSLIBS) @ATOOLSLIBS@

pkgpython_PYTHON = Sherpa.py

endif

EXTRA_DIST    = @SVNINCLUDE@ SVN_Info.C $(EXAMPLES)
