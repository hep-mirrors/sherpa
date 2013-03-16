noinst_HEADERS = \
	Sherpa.i \
	Cluster_Amplitude.i \
	Flavour.i \
	Vec4.i \
	Particle.i \
	Blob.i \
	Blob_List.i \
	Exception.i \
	Message.i

if ENABLE_PYEXT

Sherpa_wrap.cxx: $(noinst_HEADERS)
	@$(SWIG) -c++ -python -I$(top_srcdir) $<

python_LTLIBRARIES = _Sherpa.la

nodist__Sherpa_la_SOURCES = Sherpa_wrap.cxx
_Sherpa_la_CPPFLAGS= @PYTHON_CPPFLAGS@
_Sherpa_la_LDFLAGS = -module
_Sherpa_la_LIBADD  = -L@SHERPABUILDDIR@/Main -lSherpaMain \
		$(HZTOOLLIBS) $(PGSLIBS) @ATOOLSLIBS@

nodist_python_PYTHON = Sherpa.py

endif

EXTRA_DIST    = @SVNINCLUDE@
