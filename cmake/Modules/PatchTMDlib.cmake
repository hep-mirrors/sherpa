# ======================================================================
# Patches applied to the TMDlib source before configuration.
# ======================================================================

if(NOT DEFINED TMDLIB_SOURCE_DIR)
  message(FATAL_ERROR "TMDLIB_SOURCE_DIR not defined")
endif()

if(NOT DEFINED TMDLIB_STYLE_FILE)
  message(FATAL_ERROR "TMDLIB_STYLE_FILE not defined")
endif()


# ----------------------------------------------------------------------
# Patch the LaTeX style file
#
# Work around the optional sectsty dependency, which may not be
# available on the system used to build Sherpa/TMDlib.
# ----------------------------------------------------------------------

file(READ "${TMDLIB_STYLE_FILE}" _content)

string(REPLACE
  "\\usepackage{sectsty}%"
  "% \\usepackage{sectsty}%"
  _content
  "${_content}"
)

string(REPLACE
  "\\allsectionsfont{\\sffamily}"
  "% \\allsectionsfont{\\sffamily}"
  _content
  "${_content}"
)

file(WRITE "${TMDLIB_STYLE_FILE}" "${_content}")


# ----------------------------------------------------------------------
# Patch TMDlib to respect --datadir instead of hardcoding
# prefix/share/tmdlib.
# ----------------------------------------------------------------------

set(TMDLIB_MAKEFILE_IN
    "${TMDLIB_SOURCE_DIR}/include/Makefile.in")

file(READ "${TMDLIB_MAKEFILE_IN}" _makefile)

string(REPLACE
  "@prefix@/share/tmdlib/"
  "@datadir@/tmdlib/"
  _makefile
  "${_makefile}"
)

file(WRITE "${TMDLIB_MAKEFILE_IN}" "${_makefile}")

message(STATUS
  "TMDlib: patched Makefile.in to use @datadir@/tmdlib/")


# ----------------------------------------------------------------------
# Patch include ordering in src/
#
# Some external dependencies may provide broad include directories.
# Such prefixes can contain headers belonging to other packages as
# well. For example, a GSL installation prefix may also contain
# LHAPDF or APFEL++ headers.
#
# Give TMDlib's own/bundled headers and the explicitly selected
# LHAPDF headers precedence over GSL and ROOT include directories.
#
# Both Makefile.am and Makefile.in are patched. The released TMDlib
# tarball already contains Makefile.in and Sherpa invokes ./configure
# directly, without regenerating the Automake files.
#
# If the known original block is found, patch it.
# If the corrected block is already present, do nothing.
# If neither form is recognized, issue a warning and leave the file
# unchanged. This allows future TMDlib versions with a different or
# already-correct build setup to proceed normally.
# ----------------------------------------------------------------------

set(TMDLIB_SRC_MAKEFILE_AM
    "${TMDLIB_SOURCE_DIR}/src/Makefile.am")

set(TMDLIB_SRC_MAKEFILE_IN
    "${TMDLIB_SOURCE_DIR}/src/Makefile.in")


# ----------------------------------------------------------------------
# src/Makefile.am
# ----------------------------------------------------------------------

file(READ "${TMDLIB_SRC_MAKEFILE_AM}" _makefile_am)

set(_old_am
"AM_CPPFLAGS = -I$(srcdir)/../include
AM_CPPFLAGS+= $(GSL_CFLAGS)
AM_CPPFLAGS+= $(ROOT_CFLAGS)
AM_CPPFLAGS+= $(LHAPDF_CFLAGS)
AM_CPPFLAGS+= -I$(srcdir)/yaml-cpp-yaml-cpp-0.6.0/include
AM_CPPFLAGS+= -I$(srcdir)/../apfelxx/inc
AM_CPPFLAGS+= -I$(srcdir)/")

set(_new_am
"AM_CPPFLAGS = -I$(srcdir)/../include
AM_CPPFLAGS+= -I$(srcdir)/../apfelxx/inc
AM_CPPFLAGS+= -I$(srcdir)/yaml-cpp-yaml-cpp-0.6.0/include
AM_CPPFLAGS+= -I$(srcdir)/
AM_CPPFLAGS+= $(LHAPDF_CFLAGS)
AM_CPPFLAGS+= $(GSL_CFLAGS)
AM_CPPFLAGS+= $(ROOT_CFLAGS)")

string(FIND "${_makefile_am}" "${_old_am}" _am_old_pos)
string(FIND "${_makefile_am}" "${_new_am}" _am_new_pos)

if(NOT _am_old_pos EQUAL -1)

  string(REPLACE
    "${_old_am}"
    "${_new_am}"
    _makefile_am
    "${_makefile_am}"
  )

  file(WRITE
    "${TMDLIB_SRC_MAKEFILE_AM}"
    "${_makefile_am}"
  )

  message(STATUS
    "TMDlib: patched include ordering in src/Makefile.am")

elseif(NOT _am_new_pos EQUAL -1)

  message(STATUS
    "TMDlib: include ordering in src/Makefile.am is already patched")

else()

  message(WARNING
    "TMDlib: expected include-order block not found in "
    "src/Makefile.am; leaving the file unchanged")

endif()


# ----------------------------------------------------------------------
# src/Makefile.in
# ----------------------------------------------------------------------

file(READ "${TMDLIB_SRC_MAKEFILE_IN}" _makefile_in)

set(_old_in
"AM_CPPFLAGS = -I$(srcdir)/../include $(GSL_CFLAGS) $(ROOT_CFLAGS) \\
\t$(LHAPDF_CFLAGS) -I$(srcdir)/yaml-cpp-yaml-cpp-0.6.0/include \\
\t-I$(srcdir)/../apfelxx/inc -I$(srcdir)/")

set(_new_in
"AM_CPPFLAGS = -I$(srcdir)/../include \\
\t-I$(srcdir)/../apfelxx/inc \\
\t-I$(srcdir)/yaml-cpp-yaml-cpp-0.6.0/include -I$(srcdir)/ \\
\t$(LHAPDF_CFLAGS) $(GSL_CFLAGS) $(ROOT_CFLAGS)")

string(FIND "${_makefile_in}" "${_old_in}" _in_old_pos)
string(FIND "${_makefile_in}" "${_new_in}" _in_new_pos)

if(NOT _in_old_pos EQUAL -1)

  string(REPLACE
    "${_old_in}"
    "${_new_in}"
    _makefile_in
    "${_makefile_in}"
  )

  file(WRITE
    "${TMDLIB_SRC_MAKEFILE_IN}"
    "${_makefile_in}"
  )

  message(STATUS
    "TMDlib: patched include ordering in src/Makefile.in")

elseif(NOT _in_new_pos EQUAL -1)

  message(STATUS
    "TMDlib: include ordering in src/Makefile.in is already patched")

else()

  message(WARNING
    "TMDlib: expected include-order block not found in "
    "src/Makefile.in; leaving the file unchanged")

endif()


# ----------------------------------------------------------------------
# Patch include ordering in TMDplotter/ and examples-c++/
#
# These directories contain the same problematic dependency ordering:
# GSL and ROOT are placed before the explicitly selected LHAPDF, while
# TMDlib's bundled APFEL++ headers are placed last.
#
# TMDplotter and examples-c++ use the same AM_CPPFLAGS blocks, so they
# can be handled together.
#
# examples-fortran is intentionally not changed. Its targets contain
# Fortran sources and use AM_FFLAGS for their compilation.
# ----------------------------------------------------------------------

set(_old_aux_am
"AM_CPPFLAGS = -I$(srcdir)/../include 
AM_CPPFLAGS+= $(GSL_CFLAGS)
AM_CPPFLAGS+= $(ROOT_CFLAGS)
AM_CPPFLAGS+= $(LHAPDF_CFLAGS)
AM_CPPFLAGS+= -I$(top_srcdir)/src/yaml-cpp-yaml-cpp-0.6.0/include
AM_CPPFLAGS+= -I$(top_srcdir)/src/
AM_CPPFLAGS+= -I$(srcdir)/../apfelxx/inc")

set(_new_aux_am
"AM_CPPFLAGS = -I$(srcdir)/../include 
AM_CPPFLAGS+= -I$(srcdir)/../apfelxx/inc
AM_CPPFLAGS+= -I$(top_srcdir)/src/yaml-cpp-yaml-cpp-0.6.0/include
AM_CPPFLAGS+= -I$(top_srcdir)/src/
AM_CPPFLAGS+= $(LHAPDF_CFLAGS)
AM_CPPFLAGS+= $(GSL_CFLAGS)
AM_CPPFLAGS+= $(ROOT_CFLAGS)")

set(_old_aux_in
"AM_CPPFLAGS = -I$(srcdir)/../include $(GSL_CFLAGS) $(ROOT_CFLAGS) \\
\t$(LHAPDF_CFLAGS) \\
\t-I$(top_srcdir)/src/yaml-cpp-yaml-cpp-0.6.0/include \\
\t-I$(top_srcdir)/src/ -I$(srcdir)/../apfelxx/inc")

set(_new_aux_in
"AM_CPPFLAGS = -I$(srcdir)/../include \\
\t-I$(srcdir)/../apfelxx/inc \\
\t-I$(top_srcdir)/src/yaml-cpp-yaml-cpp-0.6.0/include \\
\t-I$(top_srcdir)/src/ $(LHAPDF_CFLAGS) \\
\t$(GSL_CFLAGS) $(ROOT_CFLAGS)")


foreach(_tmdlib_dir IN ITEMS TMDplotter examples-c++)

  # --------------------------------------------------------------------
  # Makefile.am
  # --------------------------------------------------------------------

  set(_makefile_am_path
      "${TMDLIB_SOURCE_DIR}/${_tmdlib_dir}/Makefile.am")

  file(READ "${_makefile_am_path}" _makefile_am)

  string(FIND "${_makefile_am}" "${_old_aux_am}" _am_old_pos)
  string(FIND "${_makefile_am}" "${_new_aux_am}" _am_new_pos)

  if(NOT _am_old_pos EQUAL -1)

    string(REPLACE
      "${_old_aux_am}"
      "${_new_aux_am}"
      _makefile_am
      "${_makefile_am}"
    )

    file(WRITE
      "${_makefile_am_path}"
      "${_makefile_am}"
    )

    message(STATUS
      "TMDlib: patched include ordering in ${_tmdlib_dir}/Makefile.am")

  elseif(NOT _am_new_pos EQUAL -1)

    message(STATUS
      "TMDlib: include ordering in ${_tmdlib_dir}/Makefile.am is already patched")

  else()

    message(WARNING
      "TMDlib: expected include-order block not found in "
      "${_tmdlib_dir}/Makefile.am; leaving the file unchanged")

  endif()


  # --------------------------------------------------------------------
  # Makefile.in
  # --------------------------------------------------------------------

  set(_makefile_in_path
      "${TMDLIB_SOURCE_DIR}/${_tmdlib_dir}/Makefile.in")

  file(READ "${_makefile_in_path}" _makefile_in)

  string(FIND "${_makefile_in}" "${_old_aux_in}" _in_old_pos)
  string(FIND "${_makefile_in}" "${_new_aux_in}" _in_new_pos)

  if(NOT _in_old_pos EQUAL -1)

    string(REPLACE
      "${_old_aux_in}"
      "${_new_aux_in}"
      _makefile_in
      "${_makefile_in}"
    )

    file(WRITE
      "${_makefile_in_path}"
      "${_makefile_in}"
    )

    message(STATUS
      "TMDlib: patched include ordering in ${_tmdlib_dir}/Makefile.in")

  elseif(NOT _in_new_pos EQUAL -1)

    message(STATUS
      "TMDlib: include ordering in ${_tmdlib_dir}/Makefile.in is already patched")

  else()

    message(WARNING
      "TMDlib: expected include-order block not found in "
      "${_tmdlib_dir}/Makefile.in; leaving the file unchanged")

  endif()

endforeach()
