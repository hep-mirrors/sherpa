pkglib_LTLIBRARIES = libExtAmp.la

GITTAG = EXTAMP
include ../git.make

libExtAmp_la_SOURCES = \
			Git_Info.C \
			External_ME_Interface.C \
			Process.C \
			Born_Process.C \
			RS_Process.C \
			BVI_Process.C \
			CS_Dipole.C \
			CS_Dipoles.C \
			Dipole_Wrapper_Process.C

if RIVET_SUPPORT
libExtAmp_la_LIBADD   =	@CONDITIONAL_RIVETLDADD@

libExtAmp_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_RIVETCPPFLAGS@
endif

localinc_HEADERS   = \
			External_ME_Interface.H \
			Process_Group.H \
			Process.H \
			Born_Process.H \
			RS_Process.H \
			BVI_Process.H \
			CS_Dipole.H \
			CS_Dipoles.H \
			Dipole_Wrapper_Process.H

