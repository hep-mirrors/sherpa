dnl workaround for old automake on darwin

AC_DEFUN([AM_CONFIG_HEADERS], [AC_CONFIG_HEADERS($@)])

dnl set flags according to build environment

AC_DEFUN([SHERPA_SETUP_BUILDSYSTEM],
[
  case "$build_os:$build_cpu:$build_vendor" in
    *darwin*:*:*)
      echo "checking for architecture... Darwin MacOS"
      if test "x$LDFLAGS" = "x"; then
        AM_LDFLAGS="-dynamic -flat_namespace"
      fi
      SEDCOMMAND="sed -i.bak -E"
      LIB_SUFFIX=".dylib"
      AC_DEFINE([ARCH_DARWIN], "1", [Architecture identified as Darwin MacOS])
      AC_DEFINE([LD_PATH_NAME], "DYLD_LIBRARY_PATH", [ld path name set to DYLD_LIBRARY_PATH]) ;;
    *linux*:*:*)
      echo "checking for architecture...  Linux"
      if test "x$LDFLAGS" = "x"; then
        AM_LDFLAGS="-rdynamic"
      fi
      SEDCOMMAND="sed -i -r"
      LIB_SUFFIX=".so"
      AC_DEFINE([ARCH_LINUX], "1", [Architecture identified as Linux])
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
    *)
      echo "checking for architecture...  unknown"
      echo "hosts system type $build not yet supported, assuming unix behaviour."
      echo "possible failure due to unknown compiler/linker characteristics."
      echo "please inform us about build results at sherpa@projects.hepforge.org"
      echo "(will continue in 10 seconds)"
      sleep 10
      if test "x$LDFLAGS" = "x"; then
        AM_LDFLAGS="-rdynamic"
      fi
      SEDCOMMAND="sed -i -r"
      LIB_SUFFIX=".so"
      AC_DEFINE([ARCH_UNIX], "1", [Architecture identified as Unix])
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
  esac
  if test "x$LDFLAGS" = "x"; then
    AX_CHECK_LINK_FLAG([-Wl,-fatal_warnings], [LDSTRICTFLAG="-fatal_warnings"], [LDSTRICTFLAG="--fatal-warnings"])
    AX_APPEND_LINK_FLAGS([-Wl,--no-as-needed], AM_LDFLAGS, [-Wl,$LDSTRICTFLAG])
  fi

  AC_DEFINE_UNQUOTED([LIB_SUFFIX], ["$LIB_SUFFIX"], [shared library suffix])
  AC_SUBST(LIB_SUFFIX)
  AC_SUBST(AM_LDFLAGS)
  if which md5sum > /dev/null; then MD5COMMAND="md5sum | cut -d' ' -f1";
  elif which openssl > /dev/null; then MD5COMMAND="openssl md5 | cut -d' ' -f2";
  else MD5COMMAND="echo 'X'"; fi
  AC_SUBST(MD5COMMAND)
  AC_SUBST(SEDCOMMAND)

  if test "x$CXXFLAGS" == "x"; then CXXFLAGS=""; fi
])


AC_DEFUN([AS_AC_EXPAND],
[
  full_var="[$2]"
  numbers="1 2 3 4"
  for i in $numbers; do
    full_var="`eval echo $full_var`";
  done
  AC_SUBST([$1], "$full_var")
])


dnl setup all variables for substitution in Makefile.am's and some additional DEFINEs
dnl
dnl Additionally some variables are defined automatically:
dnl @bindir@  executables' directory
dnl @datadir@  specified data directory e. g. /usr/local/share
dnl @includedir@  directory where header files are being installed
dnl @libdir@  directory where libraries are being installed
dnl @prefix@  the common installation prefix, e. g. /usr/local
dnl @top_builddir@  relative path to the top-level of build tree


AC_DEFUN([SHERPA_SETUP_VARIABLES],
[
  if test "x$VERSIONING" != "x"; then
    echo "x$VERSIONING";
    pkgdatadir="\${datadir}/\${PACKAGE_TARNAME}-\${VERSIONING}";
    AC_SUBST(pkgdatadir)
    pkglibdir="\${libdir}/\${PACKAGE_TARNAME}-\${VERSIONING}";
    AC_SUBST(pkglibdir)
    pkgincludedir="\${includedir}/\${PACKAGE_TARNAME}-\${VERSIONING}";
    AC_SUBST(pkgincludedir)
  else
    pkgdatadir="\${datadir}/\${PACKAGE_TARNAME}";
    AC_SUBST(pkgdatadir)
    pkglibdir="\${libdir}/\${PACKAGE_TARNAME}";
    AC_SUBST(pkglibdir)
    pkgincludedir="\${includedir}/\${PACKAGE_TARNAME}";
    AC_SUBST(pkgincludedir)
  fi;
  
  if test "x$docdir" = "x"; then
    docdir="\${datadir}/doc/\${PACKAGE_TARNAME}";
    AC_SUBST(docdir)
  fi;

  if test "x$htmldir" = "x"; then
    htmldir="\${docdir}";
    AC_SUBST(htmldir)
  fi;

  AMEGICDIR="\${top_srcdir}/AMEGIC++"
  AMEGICBUILDDIR="\${top_builddir}/AMEGIC++"
  AMEGICLIBS="\${AMEGICBUILDDIR}/Main/libAmegic.la \
	\${AMEGICBUILDDIR}/DipoleSubtraction/libDipoleSubtraction.la \
	\${AMEGICBUILDDIR}/Amplitude/libAmplitude.la \
	\${AMEGICBUILDDIR}/Phasespace/libAmegicPSGen.la \
	\${AMEGICBUILDDIR}/String/libString.la \
	\${AMEGICBUILDDIR}/Amplitude/Zfunctions/libZfunctions.la"
  AC_SUBST(AMEGICDIR)
  AC_SUBST(AMEGICBUILDDIR)
  AC_SUBST(AMEGICLIBS)

  EXTAMPDIR="\${top_srcdir}/EXTAMP"
  EXTAMPBUILDDIR="\${top_builddir}/EXTAMP"
  EXTAMPLIBS="\${EXTAMPBUILDDIR}/libExtAmp.la"
  AC_SUBST(EXTAMPDIR)
  AC_SUBST(EXTAMPBUILDDIR)
  AC_SUBST(EXTAMPLIBS)

  AMISICDIR="\${top_srcdir}/AMISIC++"
  AMISICBUILDDIR="\${top_builddir}/AMISIC++"
  AMISICLIBS="\${AMISICBUILDDIR}/Main/libAmisic.la \
        \${AMISICBUILDDIR}/Tools/libAmisicTools.la \
        \${AMISICBUILDDIR}/Perturbative/libAmisicPerturbative.la"
  AC_SUBST(AMISICDIR)
  AC_SUBST(AMISICBUILDDIR)
  AC_SUBST(AMISICLIBS)

  AHADICDIR="\${top_srcdir}/AHADIC++"
  AHADICBUILDDIR="\${top_builddir}/AHADIC++"
  AHADICLIBS="\${AHADICBUILDDIR}/Main/libAhadicMain.la \
	\${AHADICBUILDDIR}/Tools/libAhadicTools.la \
	\${AHADICBUILDDIR}/Formation/libAhadicFormation.la \
	\${AHADICBUILDDIR}/Decays/libAhadicDecays.la"
  AC_SUBST(AHADICDIR)
  AC_SUBST(AHADICBUILDDIR)
  AC_SUBST(AHADICLIBS)
  
  ATOOLSDIR="\${top_srcdir}/ATOOLS"
  ATOOLSBUILDDIR="\${top_builddir}/ATOOLS"
  ATOOLSLIBS="\${ATOOLSBUILDDIR}/Phys/libToolsPhys.la \
	\${ATOOLSBUILDDIR}/Math/libToolsMath.la \
	\${ATOOLSBUILDDIR}/YAML/libToolsYaml.la \
	\${ATOOLSBUILDDIR}/Org/libToolsOrg.la"
  AC_SUBST(ATOOLSDIR)
  AC_SUBST(ATOOLSBUILDDIR)
  AC_SUBST(ATOOLSLIBS)
  
  BEAMDIR="\${top_srcdir}/BEAM"
  BEAMBUILDDIR="\${top_builddir}/BEAM"
  BEAMLIBS="\${BEAMBUILDDIR}/Main/libBeamMain.la \
  	\${BEAMBUILDDIR}/Spectra/libBeamSpectra.la"
  AC_SUBST(BEAMDIR)
  AC_SUBST(BEAMBUILDDIR)
  AC_SUBST(BEAMLIBS)

  METOOLSDIR="\${top_srcdir}/METOOLS"
  METOOLSBUILDDIR="\${top_builddir}/METOOLS"
  METOOLSLIBS="\${METOOLSBUILDDIR}/Explicit/libMEToolsExplicit.la \
	\${METOOLSBUILDDIR}/Currents/libMEToolsCurrents.la \
	\${METOOLSBUILDDIR}/Vertices/libMEToolsVertices.la \
	\${METOOLSBUILDDIR}/Colors/libMEToolsColors.la \
	\${METOOLSBUILDDIR}/SpinCorrelations/libMEToolsSpinCorrelations.la \
	\${METOOLSBUILDDIR}/Loops/libMEToolsLoops.la \
	\${METOOLSBUILDDIR}/Main/libMEToolsMain.la"
  AC_SUBST(METOOLSDIR)
  AC_SUBST(METOOLSBUILDDIR)
  AC_SUBST(METOOLSLIBS)
  
  EXTRAXSDIR="\${top_srcdir}/EXTRA_XS"
  EXTRAXSBUILDDIR="\${top_builddir}/EXTRA_XS"
  EXTRAXSLIBS="\${EXTRAXSBUILDDIR}/Main/libExtraXS.la \
	\${EXTRAXSBUILDDIR}/Two2Two/libExtraXS2_2.la \
	\${EXTRAXSBUILDDIR}/One2Two/libExtraXS1_2.la \
	\${EXTRAXSBUILDDIR}/One2Three/libExtraXS1_3.la \
	\${EXTRAXSBUILDDIR}/NLO/libExtraXSNLO.la \
	\${EXTRAXSBUILDDIR}/Special/libExtraXSSpecial.la"
  AC_SUBST(EXTRAXSDIR)
  AC_SUBST(EXTRAXSBUILDDIR)
  AC_SUBST(EXTRAXSLIBS)
  
  MCATNLODIR="\${top_srcdir}/MCATNLO"
  MCATNLOBUILDDIR="\${top_builddir}/MCATNLO"
  MCATNLOLIBS="\${MCATNLOBUILDDIR}/Main/libMCatNLOMain.la \
	\${MCATNLOBUILDDIR}/Calculators/libMCatNLOCalculators.la \
	\${MCATNLOBUILDDIR}/Showers/libMCatNLOShowers.la \
	\${MCATNLOBUILDDIR}/Tools/libMCatNLOTools.la"
  AC_SUBST(MCATNLODIR)
  AC_SUBST(MCATNLOBUILDDIR)
  AC_SUBST(MCATNLOLIBS)
  
  CSSDIR="\${top_srcdir}/CSSHOWER++"
  CSSBUILDDIR="\${top_builddir}/CSSHOWER++"
  CSSLIBS="\${CSSBUILDDIR}/Main/libCSMain.la \
	\${CSSBUILDDIR}/Calculators/libCSCalculators.la \
	\${CSSBUILDDIR}/Showers/libCSShowers.la \
	\${CSSBUILDDIR}/Tools/libCSTools.la"
  AC_SUBST(CSSDIR)
  AC_SUBST(CSSBUILDDIR)
  AC_SUBST(CSSLIBS)
  
  DIREDIR="\${top_srcdir}/DIRE"
  DIREBUILDDIR="\${top_builddir}/DIRE"
  DIRELIBS="\${DIREBUILDDIR}/Tools/libDireTools.la \
	\${DIREBUILDDIR}/Shower/libDireShower.la \
	\${DIREBUILDDIR}/Gauge/libDireGauge.la \
	\${DIREBUILDDIR}/Lorentz/libDireLorentz.la \
	\${DIREBUILDDIR}/Main/libDireMain.la"
  AC_SUBST(DIREDIR)
  AC_SUBST(DIREBUILDDIR)
  AC_SUBST(DIRELIBS)
  
  DIMDIR="\${top_srcdir}/DIM"
  DIMBUILDDIR="\${top_builddir}/DIM"
  DIMLIBS="\${DIMBUILDDIR}/Tools/libDIMTools.la \
	\${DIMBUILDDIR}/Shower/libDIMShower.la \
	\${DIMBUILDDIR}/Gauge/libDIMGauge.la \
	\${DIMBUILDDIR}/Lorentz/libDIMLorentz.la \
	\${DIMBUILDDIR}/Main/libDIMMain.la"
  AC_SUBST(DIMDIR)
  AC_SUBST(DIMBUILDDIR)
  AC_SUBST(DIMLIBS)
  
  COMIXDIR="\${top_srcdir}/COMIX"
  COMIXBUILDDIR="\${top_builddir}/COMIX"
  COMIXLIBS="\${COMIXBUILDDIR}/Amplitude/libComixAmplitude.la \
	\${COMIXBUILDDIR}/Phasespace/libComixPhasespace.la \
	\${COMIXBUILDDIR}/Main/libComix.la"
  AC_SUBST(COMIXDIR)
  AC_SUBST(COMIXBUILDDIR)
  AC_SUBST(COMIXLIBS)
  
  HADRONSDIR="\${top_srcdir}/HADRONS++"
  HADRONSBUILDDIR="\${top_builddir}/HADRONS++"
  HADRONSLIBS="\${HADRONSBUILDDIR}/Main/libHadronsMain.la \
	\${HADRONSBUILDDIR}/ME_Library/libHadronsMEs.la \
	\${HADRONSBUILDDIR}/Current_Library/libHadronsCurrents.la \
	\${HADRONSBUILDDIR}/PS_Library/libHadronsPSs.la"
  AC_SUBST(HADRONSDIR)
  AC_SUBST(HADRONSBUILDDIR)
  AC_SUBST(HADRONSLIBS)
  
  PHOTONSDIR="\${top_srcdir}/PHOTONS++"
  PHOTONSBUILDDIR="\${top_builddir}/PHOTONS++"
  PHOTONSLIBS="\${PHOTONSBUILDDIR}/Main/libPhotonsMain.la \
	\${PHOTONSBUILDDIR}/Tools/libPhotonsTools.la \
	\${PHOTONSBUILDDIR}/PhaseSpace/libPhotonsPhaseSpace.la \
	\${PHOTONSBUILDDIR}/MEs/libPhotonsMEs.la"
  AC_SUBST(PHOTONSDIR)
  AC_SUBST(PHOTONSBUILDDIR)
  AC_SUBST(PHOTONSLIBS)
  
  MODELDIR="\${top_srcdir}/MODEL"
  MODELBUILDDIR="\${top_builddir}/MODEL"
  MODELMAINLIB="\${MODELBUILDDIR}/Main/libModelMain.la"
  MODELLIBS="\${MODELMAINLIB} \
	\${MODELBUILDDIR}/UFO/libModelUFO.la"
  AC_SUBST(MODELDIR)
  AC_SUBST(MODELBUILDDIR)
  AC_SUBST(MODELMAINLIB)
  AC_SUBST(MODELLIBS)
  
  PDFDIR="\${top_srcdir}/PDF"
  PDFBUILDDIR="\${top_builddir}/PDF"
  PDFINCS="-I\${PDFDIR}/Main"
  PDFLIBS="\${PDFBUILDDIR}/Main/libPDF.la"
  AC_SUBST(PDFDIR)
  AC_SUBST(PDFBUILDDIR)
  AC_SUBST(PDFLIBS)
  
  REMNANTSDIR="\${top_srcdir}/REMNANTS"
  REMNANTSBUILDDIR="\${top_builddir}/REMNANTS"
  REMNANTSLIBS="\${REMNANTSBUILDDIR}/Tools/libRemnantsTools.la \
        \${REMNANTSBUILDDIR}/Main/libRemnants.la"
  AC_SUBST(REMNANTSDIR)
  AC_SUBST(REMNANTSBUILDDIR)
  AC_SUBST(REMNANTSLIBS)

  RECONNECTIONSDIR="\${top_srcdir}/RECONNECTIONS"
  RECONNECTIONSBUILDDIR="\${top_builddir}/RECONNECTIONS"
  RECONNECTIONSLIBS="\${RECONNECTIONSBUILDDIR}/Main/libReconnections.la"
  AC_SUBST(RECONNECTIONSDIR)
  AC_SUBST(RECONNECTIONSBUILDDIR)
  AC_SUBST(RECONNECTIONSLIBS)

  PHASICDIR="\${top_srcdir}/PHASIC++"
  PHASICBUILDDIR="\${top_builddir}/PHASIC++"
  PHASICLIBS="\${PHASICBUILDDIR}/Main/libPhasicMain.la \
	\${PHASICBUILDDIR}/Channels/libPhasicChannels.la \
	\${PHASICBUILDDIR}/Process/libPhasicProcess.la \
	\${PHASICBUILDDIR}/Selectors/libPhasicSelectors.la \
	\${PHASICBUILDDIR}/Scales/libPhasicScales.la \
	\${PHASICBUILDDIR}/Enhance/libPhasicEnhance.la \
	\${PHASICBUILDDIR}/Decays/libPhasicDecays.la"
  AC_SUBST(PHASICDIR)
  AC_SUBST(PHASICBUILDDIR)
  AC_SUBST(PHASICLIBS)
  
  SHRIMPSDIR="\${top_srcdir}/SHRiMPS"
  SHRIMPSBUILDDIR="\${top_builddir}/SHRiMPS"
  SHRIMPSLIBS="\${SHRIMPSBUILDDIR}/Main/libShrimpsMain.la \
	\${SHRIMPSBUILDDIR}/Event_Generation/libShrimpsEvents.la \
	\${SHRIMPSBUILDDIR}/Beam_Remnants/libShrimpsBeamRemnants.la \
	\${SHRIMPSBUILDDIR}/Cross_Sections/libShrimpsXsecs.la \
	\${SHRIMPSBUILDDIR}/Eikonals/libShrimpsEikonals.la \
	\${SHRIMPSBUILDDIR}/Tools/libShrimpsTools.la"
  AC_SUBST(SHRIMPSDIR)
  AC_SUBST(SHRIMPSBUILDDIR)
  AC_SUBST(SHRIMPSLIBS)

  SHERPADIR="\${top_srcdir}/SHERPA"
  SHERPABUILDDIR="\${top_builddir}/SHERPA"
  SHERPALIBS="\${SHERPABUILDDIR}/Initialization/libSherpaInitialization.la \
	\${SHERPABUILDDIR}/Single_Events/libSherpaSingleEvents.la \
	\${SHERPABUILDDIR}/SoftPhysics/libSherpaSoftPhysics.la \
	\${SHERPABUILDDIR}/PerturbativePhysics/libSherpaPerturbativePhysics.la \
	\${SHERPABUILDDIR}/LundTools/libLundTools.la \
	\${SHERPABUILDDIR}/Tools/libSherpaTools.la"
  AC_SUBST(SHERPADIR)
  AC_SUBST(SHERPABUILDDIR)
  AC_SUBST(SHERPALIBS)

  if test "x$prefix" = "xNONE"; then
    prefix=$ac_default_prefix
  fi
  if test "x$exec_prefix" = "xNONE"; then
    exec_prefix=$prefix
  fi

  AS_AC_EXPAND(LIBDIR, ${pkglibdir})
  AS_AC_EXPAND(PYLIBDIR, ${pythondir})
  AS_AC_EXPAND(INCLUDEDIR, ${pkgincludedir})
  AS_AC_EXPAND(BINDIR, ${bindir})
  AS_AC_EXPAND(DATADIR, ${pkgdatadir})
  AS_AC_EXPAND(SHERPAPREFIX, ${prefix})
  AS_AC_EXPAND(PYTHONLIBS, ${pythondir})

  AC_DEFINE_UNQUOTED([SHERPA_VERSION], ["`echo AC_PACKAGE_VERSION | cut -d. -f1`"], [Sherpa version])
  AC_DEFINE_UNQUOTED([SHERPA_SUBVERSION], ["`echo AC_PACKAGE_VERSION | cut -d. -f2,3`"], [Sherpa subversion])
  AC_DEFINE_UNQUOTED([SHERPA_PREFIX], "$SHERPAPREFIX", [Sherpa installation prefix])
  AC_DEFINE_UNQUOTED([SHERPA_INCLUDE_PATH], "$INCLUDEDIR", [Sherpa include directory])
  AC_DEFINE_UNQUOTED([SHERPA_LIBRARY_PATH], "$LIBDIR", [Sherpa library directory])
  AC_DEFINE_UNQUOTED([SHERPA_SHARE_PATH], "$DATADIR", [Sherpa data directory])
  AC_DEFINE_UNQUOTED([PYTHON_LIBS], "$PYTHONLIBS", [Sherpa python library directory])
  AC_DEFINE([USING__COLOUR], "1", [Using colour])

  AM_CPPFLAGS="-I\$(top_srcdir)"
  AC_SUBST(AM_CPPFLAGS)

  AM_CXXFLAGS="-g -O2"
  AC_LANG_PUSH([C++])
  AX_CHECK_COMPILE_FLAG(
    -fcx-fortran-rules,
    [AM_CXXFLAGS="${AM_CXXFLAGS} -fcx-fortran-rules"])
  AC_LANG_POP([C++])
  AC_SUBST(AM_CXXFLAGS)

  localincdir="\$(pkgincludedir)/\$(subdir)"
  AC_SUBST(localincdir)
])



dnl Conditional compiling and linking

AC_DEFUN([SHERPA_SETUP_CONFIGURE_OPTIONS],
[
  AC_ARG_ENABLE(
    versioning,
    AS_HELP_STRING([--enable-versioning],[Add version tag to executables and library/header directories, such that multiple Sherpa versions can live in the same prefix.]),
    [ AC_MSG_CHECKING(whether to enable versioning)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no);
             VERSIONING="";;
        yes) AC_MSG_RESULT(yes);
             VERSIONING="AC_PACKAGE_VERSION";;
        *)   if test "x${enableval}" != "x"; then
               AC_MSG_RESULT(yes);
               VERSIONING="${enableval}"
             fi
      esac ],
    [ AC_MSG_CHECKING(whether to enable versioning);
      AC_MSG_RESULT(no);
      VERSIONING=""; ] 
  )
  AC_SUBST([VERSIONING])

  AC_ARG_ENABLE(
    multithread,
    AS_HELP_STRING([--enable-multithread],[Enable multithreading]),
    [ AC_MSG_CHECKING(for multithreading)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); multithread=false ;;
        yes) AC_MSG_RESULT(yes); multithread=true ;;
      esac ],
    [ AC_MSG_CHECKING(for multithreading); AC_MSG_RESULT(no); multithread=false ] 
  )
  if test "$multithread" = "true" ; then
    AC_DEFINE([USING__Threading], "1", [using multithreading])
    CONDITIONAL_THREADLIBS="-lpthread"
  fi
  AC_SUBST(CONDITIONAL_THREADLIBS)
  AM_CONDITIONAL(USING__Threading, test "$multithread" = "true" )
  
  AC_ARG_ENABLE(
    analysis,
    AS_HELP_STRING([--enable-analysis],[Enable analysis]),
    [ AC_MSG_CHECKING(for analysis)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); analysis=false ;;
        yes) AC_MSG_RESULT(yes); analysis=true ;;
      esac ],
    [ AC_MSG_CHECKING(for analysis); AC_MSG_RESULT(no); analysis=false ]
  )
  AM_CONDITIONAL(USING__Analysis, test "$analysis" = "true" )

  AC_ARG_ENABLE(
    root,
    AS_HELP_STRING([--enable-root=/path/to/root],[Enable ROOT support and specify where it is installed if non-standard.]),
    [ AC_MSG_CHECKING(for ROOT installation directory)
      case "${enableval}" in
        no)  AC_MSG_RESULT(ROOT not enabled); root=false;;
        yes) if test -d "$ROOTSYS"; then
               CONDITIONAL_ROOTDIR=$ROOTSYS
               CONDITIONAL_ROOTINCS="-I$ROOTSYS/include -I$($ROOTSYS/bin/root-config --incdir)";
               CONDITIONAL_ROOTLIBS="-L$ROOTSYS/lib $($ROOTSYS/bin/root-config --glibs)"
               CONDITIONAL_ROOTFLAGS="$($ROOTSYS/bin/root-config --cflags)"
             elif test -x "`which root-config`"; then
               CONDITIONAL_ROOTDIR=`root-config --prefix`;
               CONDITIONAL_ROOTINCS=-I`root-config --incdir`;
               CONDITIONAL_ROOTLIBS=`root-config --glibs`;
               CONDITIONAL_ROOTFLAGS=`root-config --cflags`;
                if ! test -d "$CONDITIONAL_ROOTDIR"; then
                  AC_MSG_ERROR(root-config --prefix returned a path that is not available. Please check your ROOT installation and set \$ROOTSYS manually.);
                fi
             else
               AC_MSG_ERROR(\$ROOTSYS is not a valid path and root-config was not found.);
             fi;
             AC_MSG_RESULT([${CONDITIONAL_ROOTDIR}]); root=true;;
        *)   if test -d "${enableval}"; then
               CONDITIONAL_ROOTDIR="${enableval}"
               CONDITIONAL_ROOTINCS="-I${enableval}/include -I${enableval}/include/root";
               CONDITIONAL_ROOTLIBS="-L${enableval}/lib $(${enableval}/bin/root-config --glibs)";
               CONDITIONAL_ROOTFLAGS="$(${enableval}/bin/root-config --cflags)";
             else
               AC_MSG_ERROR(${enableval} is not a valid path.);
             fi;
             AC_MSG_RESULT([${CONDITIONAL_ROOTDIR}]); root=true;;
      esac ],
    [ root=false ]
  )
  if test "$root" = "true" ; then
    AC_DEFINE([USING__ROOT], "1", [using ROOT])
    fi
  AC_SUBST(CONDITIONAL_ROOTDIR)
  AC_SUBST(CONDITIONAL_ROOTINCS)
  AC_SUBST(CONDITIONAL_ROOTLIBS)
  AC_SUBST(CONDITIONAL_ROOTFLAGS)
  AM_CONDITIONAL(ROOT_SUPPORT, test "$root" = "true")
  

  AC_ARG_ENABLE(
    hepmc2,
    AS_HELP_STRING([--enable-hepmc2=/path/to/hepmc],[Enable HepMC (version 2.x) support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for HepMC2 installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(HepMC2 not enabled); hepmc2=false ;;
        yes)  if test -d "$HEPMC2DIR"; then
                CONDITIONAL_HEPMC2DIR="$HEPMC2DIR"
                CONDITIONAL_HEPMC2INCS="-I$HEPMC2DIR/include"
                CONDITIONAL_HEPMC2LIBS="-L$HEPMC2DIR/lib -Wl,-rpath -Wl,$HEPMC2DIR/lib -L$HEPMC2DIR/lib64 -Wl,-rpath -Wl,$HEPMC2DIR/lib64 -lHepMC";
              else
                AC_MSG_ERROR(\$HEPMC2DIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_HEPMC2DIR}]); hepmc2=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_HEPMC2DIR="${enableval}"
                CONDITIONAL_HEPMC2INCS="-I${enableval}/include"
                CONDITIONAL_HEPMC2LIBS="-L${enableval}/lib -Wl,-rpath -Wl,${enableval}/lib -L${enableval}/lib64 -Wl,-rpath -Wl,${enableval}/lib64 -lHepMC";
              else
                AC_MSG_ERROR(${enableval} is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_HEPMC2DIR}]); hepmc2=true;;
      esac
      if test -f "$CONDITIONAL_HEPMC2DIR/include/HepMC/IO_GenEvent.h"; then
        hepmciogenevent=true;
      fi;
      if test -f "$CONDITIONAL_HEPMC2DIR/include/HepMC/HepMCDefs.h"; then
        hepmcdefs=true;
      fi;
      if test -f "$CONDITIONAL_HEPMC2DIR/include/HepMC/Units.h"; then
        hepmcunits=true;
      fi;
      ],
    [ hepmc2=false ]
  )
  if test "$hepmc2" = "true" ; then
    AC_DEFINE([USING__HEPMC2], "1", [Using HEPMC2])
    if test "$hepmciogenevent" = "true"; then
      AC_DEFINE([USING__HEPMC2__IOGENEVENT], "1", [HepMC::IO_GenEvent available])
    fi
    if test "$hepmcunits" = "true"; then
      AC_DEFINE([USING__HEPMC2__UNITS], "1", [HepMC::Units available])
    fi
    if test "$hepmcdefs" = "true"; then
      AC_DEFINE([USING__HEPMC2__DEFS], "1", [HepMCDefs.h available])
    fi
  fi
  AC_SUBST(CONDITIONAL_HEPMC2DIR)
  AC_SUBST(CONDITIONAL_HEPMC2INCS)
  AC_SUBST(CONDITIONAL_HEPMC2LIBS)
  AM_CONDITIONAL(HEPMC2_SUPPORT, test "$hepmc2" = "true")

  AC_ARG_ENABLE(
    hepmc3root,
    AS_HELP_STRING([--enable-hepmc3root],[Enable HepMC (version 3.1+) ROOT support]),
    [ 
    case "${enableval}" in
        no)  AC_MSG_RESULT(HepMC3 ROOT support not enabled); hepmc3root=false ;;
        yes) AC_MSG_RESULT(HepMC3 ROOT support enabled); hepmc3root=true ;;
    
          esac
    ],
    [ hepmc3root=${root} ]
  )


  AC_ARG_ENABLE(
    hepmc3,
    AS_HELP_STRING([--enable-hepmc3=/path/to/hepmc],[Enable HepMC (version 3.x) support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for HepMC3 installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(HepMC3 not enabled); hepmc3=false  ;;
        yes) if test -x "`which HepMC3-config`"; then
               CONDITIONAL_HEPMC3DIR=`HepMC3-config --prefix`;
             fi;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_HEPMC3DIR=${enableval};
            fi;;
      esac;
      if test -n "$CONDITIONAL_HEPMC3DIR"; then
        if test -x "$CONDITIONAL_HEPMC3DIR/bin/HepMC3-config"; then      
          AC_MSG_RESULT([${CONDITIONAL_HEPMC3DIR}]); hepmc3=true
          CONDITIONAL_HEPMC3INCS="$($CONDITIONAL_HEPMC3DIR/bin/HepMC3-config --cppflags)";
          CONDITIONAL_HEPMC3LIBS="$($CONDITIONAL_HEPMC3DIR/bin/HepMC3-config --libs)";
          if test "$hepmc3root" = "true" ; then
            CONDITIONAL_HEPMC3INCS="$($CONDITIONAL_HEPMC3DIR/bin/HepMC3-config --cppflags --rootIO) ${CONDITIONAL_ROOTINCS}";
            CONDITIONAL_HEPMC3LIBS="$($CONDITIONAL_HEPMC3DIR/bin/HepMC3-config --libs  --rootIO) ${CONDITIONAL_ROOTLIBS}";
            if ! test -f "$($CONDITIONAL_HEPMC3DIR/bin/HepMC3-config --includedir)/HepMC3/WriterRoot.h" -a -f "$($CONDITIONAL_HEPMC3DIR/bin/HepMC3-config --includedir)/HepMC3/WriterRootTree.h"; then
               AC_MSG_ERROR(HepMC3 installation does not contain ROOT support.);
            fi;
          fi;
        else
          AC_MSG_ERROR(Unable to use HepMC3 from specified path);
        fi;
      fi;
    ],
    [ hepmc3=false ]
  )
  if test "$hepmc3" = "true" ; then
    AC_DEFINE([USING__HEPMC3], "1", [Using HEPMC3])
  fi
  if test "$hepmc3root" = "true" ; then
    AC_DEFINE([USING__HEPMC3__ROOT], "1", [HepMC3 with ROOT support])
  fi
  AC_SUBST(CONDITIONAL_HEPMC3DIR)
  AC_SUBST(CONDITIONAL_HEPMC3INCS)
  AC_SUBST(CONDITIONAL_HEPMC3LIBS)
  AM_CONDITIONAL(HEPMC3_SUPPORT, test "$hepmc3" = "true")


  AC_ARG_ENABLE(
    rivet,
    AS_HELP_STRING([--enable-rivet=/path/to/rivet],[Enable Rivet support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for Rivet installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(Rivet not enabled); rivet2=false; rivet3=false ;;
        yes) if test -x "`which rivet-config`"; then
               CONDITIONAL_RIVETDIR=`rivet-config --prefix`;
             fi;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_RIVETDIR=${enableval};
            fi;;
      esac;
      if test -n "$CONDITIONAL_RIVETDIR"; then
        if test -x "$CONDITIONAL_RIVETDIR/bin/rivet-config"; then
          CONDITIONAL_RIVETLDADD="$($CONDITIONAL_RIVETDIR/bin/rivet-config --ldflags) $($CONDITIONAL_RIVETDIR/bin/rivet-config --ldadd)";
          CONDITIONAL_RIVETCPPFLAGS="$($CONDITIONAL_RIVETDIR/bin/rivet-config --cppflags)";
          AC_MSG_RESULT([${CONDITIONAL_RIVETDIR}]);
          rivetversion="$($CONDITIONAL_RIVETDIR/bin/rivet-config --version)"
          AC_MSG_CHECKING(for Rivet version)
          AX_COMPARE_VERSION([${rivetversion}],[ge],[3.1.1],[ rivet3=true; AC_MSG_RESULT(Rivet 3) ], [
            AX_COMPARE_VERSION([${rivetversion}],[ge],[3.0.0],[ AC_MSG_ERROR(Rivet version 3.0.0-3.1.0 not supported -- please use 3.1.1 or above.) ], [
              AX_COMPARE_VERSION([${rivetversion}],[ge],[2.0.0],[ rivet2=true; AC_MSG_RESULT(Rivet 2) ], [
                AC_MSG_ERROR(Rivet version <2.0 found, not supported.)
              ])
            ])
          ])
        else
          AC_MSG_ERROR(Unable to use Rivet from specified path.);
        fi;
      fi;
    ],
    [ rivet=false ]
  )
  AC_SUBST(CONDITIONAL_RIVETLDADD)
  AC_SUBST(CONDITIONAL_RIVETCPPFLAGS)
  if test "$rivet2" = "true" ; then
    AC_DEFINE([USING__RIVET2], "1", [using Rivet2])
  fi
  if test "$rivet3" = "true" ; then
    AC_DEFINE([USING__RIVET3], "1", [using Rivet3])
  fi
  AM_CONDITIONAL(RIVET_SUPPORT, test "$rivet2" = "true" -o "$rivet3" = "true")

  AC_ARG_ENABLE([manual],
    AS_HELP_STRING([--enable-manual], [Enable the manual]),
      [ AC_MSG_CHECKING(whether the manual dependencies are installed);
      case "${enableval}" in
        no)  AC_MSG_RESULT(Manual not enabled); manual=false ;;
        yes)  if ! command -v python3 &>/dev/null; then
                 AC_MSG_ERROR(python3 not installed.);
              fi;

              if ! command -v sphinx-build &>/dev/null; then
                  AC_MSG_ERROR(sphinx not installed.);
              fi;

              if ! command -v makeinfo &>/dev/null; then
                  AC_MSG_ERROR(makeinfo not installed.);
              fi;

              if ! command -v pdflatex --version &>/dev/null; then
                  AC_MSG_ERROR(pdflatex not installed.);
              fi;

              # Check sphinx version number: > __MA_REQ.x.x
              __VERSION=$(sphinx-build --version | cut -d ' ' -f 2)
              __MA_VERS=$(echo $VERSION | cut -d '.' -f 1)
              __MA_REQ=2

              if !  [[ "$__MA_VERS" -ge "$__MA_REQ" ]] ; then
                  AC_MSG_ERROR("Sphinx version >= $__MA_REQ.x required. You have version: $__VERSION");
              fi;

              if ! [ python3 -c  'import sphinxcontrib.bibtex' &>/dev/null ]; then
                 AC_MSG_ERROR(sphinxcontrib-bibtex not installed.);
              fi;

              AC_MSG_RESULT([yes]); manual=true;;
           *) ;;
      esac
      ],
      [ manual=false ])

  AM_CONDITIONAL(Manual_ENABLED, test "$manual" = "true")

  AC_ARG_ENABLE(
    blackhat,
    AS_HELP_STRING([--enable-blackhat=/path/to/blackhat],[Enable BLACKHAT.]),
    [ AC_MSG_CHECKING(for BLACKHAT installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(BLACKHAT not enabled); blackhat=false ;;
        yes)  if test -x "$BLACKHATDIR/bin/blackhat-config"; then
                CONDITIONAL_BLACKHATDIR="$BLACKHATDIR"
                CONDITIONAL_BLACKHATINCS="-I$($CONDITIONAL_BLACKHATDIR/bin/blackhat-config --include)";
                CONDITIONAL_BLACKHATLIBS="$($CONDITIONAL_BLACKHATDIR/bin/blackhat-config --libs)"
              else
                AC_MSG_ERROR(\$BLACKHATDIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_BLACKHATDIR}]); blackhat=true;;
        *)    if test -x "${enableval}/bin/blackhat-config"; then
                CONDITIONAL_BLACKHATDIR="${enableval}"
                CONDITIONAL_BLACKHATINCS="-I$($CONDITIONAL_BLACKHATDIR/bin/blackhat-config --include)";
                CONDITIONAL_BLACKHATLIBS="$($CONDITIONAL_BLACKHATDIR/bin/blackhat-config --libs)"
              else
                AC_MSG_ERROR(${enableval} is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_BLACKHATDIR}]); blackhat=true;;
      esac
      ],
    [ blackhat=false ]
  )
  if test "$blackhat" = "true" ; then
    AC_DEFINE_UNQUOTED([BLACKHAT_PATH], "$CONDITIONAL_BLACKHATDIR", [BlackHat directory])
    AC_DEFINE([USING__BLACKHAT], "1", [Using BLACKHAT])
  fi
  AC_SUBST(CONDITIONAL_BLACKHATDIR)
  AC_SUBST(CONDITIONAL_BLACKHATINCS)
  AC_SUBST(CONDITIONAL_BLACKHATLIBS)
  AM_CONDITIONAL(BLACKHAT_SUPPORT, test "$blackhat" = "true")

  AC_ARG_ENABLE(
    openloops,
    AS_HELP_STRING([--enable-openloops=/path/to/openloops],[Enable OpenLoops.]),
    [ AC_MSG_CHECKING(for OpenLoops installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(OpenLoops not enabled); openloops=false ;;
        *)   OPENLOOPS_PREFIX="$(echo ${enableval} | sed -e 's/\/$//g')"
             openloops=true;
             if test -d "${OPENLOOPS_PREFIX}"; then
                AC_MSG_RESULT([${OPENLOOPS_PREFIX}]);
             else
                AC_MSG_WARN(${OPENLOOPS_PREFIX} is not a valid path.);
             fi;;
      esac
      ],
    [ openloops=false ]
  )
  if test "$openloops" = "true" ; then
    AC_DEFINE_UNQUOTED([OPENLOOPS_PREFIX], "$OPENLOOPS_PREFIX", [Openloops installation prefix])
  fi
  AM_CONDITIONAL(OPENLOOPS_SUPPORT, test "$openloops" = "true")

  
  AC_ARG_ENABLE(
    recola,
    AC_HELP_STRING([--enable-recola=/path/to/recola], [Enable Recola.]),
    [ AC_MSG_CHECKING(for Recola installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(Recola not enabled); recola=false ;;
        *)   RECOLA_PREFIX="$(echo ${enableval} | sed -e 's/\/$//g')"
             recola=true;
             if test -d "${RECOLA_PREFIX}"; then
                AC_MSG_RESULT([${RECOLA_PREFIX}]);
		CONDITIONAL_RECOLAINCS="-I$RECOLA_PREFIX/include";
             else
                AC_MSG_WARN(${RECOLA_PREFIX} is not a valid path.);
             fi;;
      esac
      ],
    [ recola=false ]
  )
  if test "$recola" = "true" ; then
    AC_DEFINE_UNQUOTED([RECOLA_PREFIX], "$RECOLA_PREFIX", [Recola installation prefix])
  fi
  AC_SUBST(CONDITIONAL_RECOLAINCS)
  AM_CONDITIONAL(RECOLA_SUPPORT, test "$recola" = "true")


  AC_ARG_ENABLE(
    gosam,
    AS_HELP_STRING([--enable-gosam=/path/to/gosam],[Enable GoSam.]),
    [ AC_MSG_CHECKING(for GoSam installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(GoSam not enabled); gosam=false ;;
        *)   GOSAM_PREFIX="$(echo ${enableval} | sed -e 's/\/$//g')"
             gosam=true;
             if test -d "${GOSAM_PREFIX}"; then
                AC_MSG_RESULT([${GOSAM_PREFIX}]);
             else
                AC_MSG_WARN(${GOSAM_PREFIX} is not a valid path.);
             fi;;
      esac
      ],
    [ gosam=false ]
  )
  if test "$gosam" = "true" ; then
    AC_DEFINE_UNQUOTED([GOSAM_PREFIX], "$GOSAM_PREFIX", [GoSam installation prefix])
  fi
  AM_CONDITIONAL(GOSAM_SUPPORT, test "$gosam" = "true")

  AC_ARG_ENABLE(
    madloop,
    AS_HELP_STRING([--enable-madloop=/path/to/madloop],[Enable Madloop.]),
    [ AC_MSG_CHECKING(for Madloop installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(Madloop not enabled); madloop=false ;;
        *)   MADLOOP_PREFIX="$(echo ${enableval} | sed -e 's/\/$//g')"
             madloop=true;
             if test -d "${MADLOOP_PREFIX}"; then
                AC_MSG_RESULT([${MADLOOP_PREFIX}]);
             else
                AC_MSG_WARN(${MADLOOP_PREFIX} is not a valid path.);
             fi;;
      esac
      ],
    [ madloop=false ]
  )
  if test "$madloop" = "true" ; then
    AC_DEFINE_UNQUOTED([MADLOOP_PREFIX], "$MADLOOP_PREFIX", [Madloop installation prefix])
  fi
  AM_CONDITIONAL(MADLOOP_SUPPORT, test "$madloop" = "true")
  
  AC_ARG_ENABLE(
    mcfm,
    AS_HELP_STRING([--enable-mcfm=/path/to/mcfm],[Enable MCFM.]),
    [ AC_MSG_CHECKING(for MCFM installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(MCFM not enabled); mcfm=false ;;
        yes)  if test -d "$MCFMDIR"; then
                CONDITIONAL_MCFMDIR="$MCFMDIR"
                CONDITIONAL_MCFMLIBS="-Wl,-rpath -Wl,$CONDITIONAL_MCFMDIR/lib -L$CONDITIONAL_MCFMDIR/lib -lMCFM"
                CONDITIONAL_MCFMINCS="-I$CONDITIONAL_MCFMDIR/include"
              else
                AC_MSG_ERROR(\$MCFMDIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_MCFMDIR}]); mcfm=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_MCFMDIR="${enableval}"
                CONDITIONAL_MCFMLIBS="-Wl,-rpath -Wl,$CONDITIONAL_MCFMDIR/lib -L$CONDITIONAL_MCFMDIR/lib -lMCFM"
                CONDITIONAL_MCFMINCS="-I$CONDITIONAL_MCFMDIR/include"
              else
                AC_MSG_ERROR(${enableval} is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_MCFMDIR}]); mcfm=true;;
      esac
      ],
    [ mcfm=false ]
  )
  if test "$mcfm" = "true" ; then
    AC_DEFINE([USING__MCFM], "1", [Using MCFM])
    AC_DEFINE_UNQUOTED([MCFM_PATH], "$CONDITIONAL_MCFMDIR", [MCFM directory])
  fi
  AC_SUBST(CONDITIONAL_MCFMDIR)
  AC_SUBST(CONDITIONAL_MCFMLIBS)
  AC_SUBST(CONDITIONAL_MCFMINCS)
  AM_CONDITIONAL(MCFM_SUPPORT, test "$mcfm" = "true")

  AC_ARG_ENABLE(
    lhole,
    AS_HELP_STRING([--enable-lhole],[Enable Les Houches One-Loop Generator interface.]),
    [ AC_MSG_CHECKING(for LHOLE)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); lhole=false ;;
        yes) AC_MSG_RESULT(yes); lhole=true ;;
      esac ],
    [ AC_MSG_CHECKING(for LHOLE); AC_MSG_RESULT(no); lhole=false ]
  )
  AM_CONDITIONAL(USING__LHOLE, test "$lhole" = "true" )

  AC_ARG_WITH(
    lhapdf,
    AS_HELP_STRING([--with-lhapdf=@<:@ARG@:>@],[use LHAPDF from @<:@ARG@:>@.]),
    [ AC_MSG_CHECKING(for LHAPDF installation directory);
      if test "$withval" = "install"; then
        CONDITIONAL_LHAPDFDIR=$ac_default_prefix;
        lhapdf_version="6.4.0";
        test "x$prefix" != xNONE && CONDITIONAL_LHAPDFDIR=$prefix;
        if ! test -f ${CONDITIONAL_LHAPDFDIR}/bin/lhapdf-config; then
          wget --no-check-certificate https://lhapdf.hepforge.org/downloads/LHAPDF-${lhapdf_version}.tar.gz
          tar xzf LHAPDF-${lhapdf_version}.tar.gz;
          cd LHAPDF-${lhapdf_version};
          ./configure --prefix=${CONDITIONAL_LHAPDFDIR} --disable-python || exit;
          make || exit; make install || exit;
          cd ${CONDITIONAL_LHAPDFDIR}/share/LHAPDF;
          wget --no-check-certificate https://www.hep.ucl.ac.uk/pdf4lhc/pdf4lhc21grids/PDF4LHC21_40_pdfas.tgz;
          tar xzf PDF4LHC21_40_pdfas.tgz;
          echo "93300 PDF4LHC21_40_pdfas 1" >> pdfsets.index
          cd -;
          cd ..;
          rm -rf LHAPDF-${lhapdf_version}.tar.gz LHAPDF-${lhapdf_version} ${CONDITIONAL_LHAPDFDIR}/share/LHAPDF/PDF4LHC21_40_pdfas.tgz;
          echo "Successfully installed lhapdf into ${CONDITIONAL_LHAPDFDIR}."
        fi
      else
        CONDITIONAL_LHAPDFDIR="$withval"
      fi
    ],
    [ CONDITIONAL_LHAPDFDIR=/usr; ]
  )
  if test -x "$CONDITIONAL_LHAPDFDIR/bin/lhapdf-config"; then
    CONDITIONAL_LHAPDFLIBS="$($CONDITIONAL_LHAPDFDIR/bin/lhapdf-config --ldflags)";
    CONDITIONAL_LHAPDFINCS="$($CONDITIONAL_LHAPDFDIR/bin/lhapdf-config --cppflags)";
    AC_MSG_RESULT([${CONDITIONAL_LHAPDFDIR}]); lhapdf=true;
  else
    AC_MSG_ERROR(Did not find required dependency LHAPDF in ${CONDITIONAL_LHAPDFDIR}. Please specify its installation prefix using '--with-lhapdf=/path' or enable its automatic installation using '--with-lhapdf=install');
  fi;
  AC_SUBST(CONDITIONAL_LHAPDFDIR)
  AC_SUBST(CONDITIONAL_LHAPDFLIBS)
  AC_SUBST(CONDITIONAL_LHAPDFINCS)
  AC_DEFINE_UNQUOTED([LHAPDF_PATH], "$CONDITIONAL_LHAPDFDIR", [LHAPDF directory])

  AC_ARG_ENABLE(
    hztool,
    AS_HELP_STRING([--enable-hztool=/path/to/hztool],[Enable hztool for analysis.]),
    [ AC_MSG_CHECKING(for hztool installation directory);
      case "${enableval}" in
        no) AC_MSG_RESULT(hztool not enabled); hztool=false;;
        *)  if test -d "${enableval}"; then
              if test -f "${enableval}/lib/libhztool.so"; then
                CONDITIONAL_HZTOOLLIBS="-L${enableval}/lib -lhztool";
                CONDITIONAL_HZTOOLINCS="-I${enableval}/include/hztool";
                CONDITIONAL_HZTOOLDIR="${enableval}";
                hztool=true;
                AC_MSG_RESULT(${enableval});
              else
              if test -f "${enableval}/lib64/libhztool.so"; then
                CONDITIONAL_HZTOOLLIBS="-L${enableval}/lib64 -lhztool";
                CONDITIONAL_HZTOOLINCS="-I${enableval}/include/hztool";
                CONDITIONAL_HZTOOLDIR="${enableval}";
                hztool=true;
                AC_MSG_RESULT(${enableval});
               else
                AC_MSG_ERROR(Did not find '${enableval}/lib/libhztool.so and ${enableval}/lib64/libhztool.so'.); 
              fi;
              fi;
            else
              AC_MSG_ERROR(Did not find hztool directory '${enableval}'.);
            fi;
      esac;
    ],
    [ hztool=false ]
  )
  if test "$hztool" = "true" ; then
    AC_DEFINE([USING__HZTOOL], "1", [hztool found])
  fi
  AC_SUBST(CONDITIONAL_HZTOOLDIR)
  AC_SUBST(CONDITIONAL_HZTOOLINCS)
  AC_SUBST(CONDITIONAL_HZTOOLLIBS)
  AM_CONDITIONAL(HZTOOL_SUPPORT, test "$hztool" = "true")

  AC_ARG_ENABLE(
    cernlib,
    AS_HELP_STRING([--enable-cernlib=/path/to/cernlib],[Enable cernlib.]),
    [ AC_MSG_CHECKING(for cernlib installation directory);
      case "${enableval}" in
        no) AC_MSG_RESULT(cernlib not enabled); cernlib=false;;
        *)  if test -d "${enableval}"; then
              if test -f "${enableval}/lib/libkernlib_noshift.a"; then
                if test -f "${enableval}/lib/libkernlib_noshift.so"; then
                  CONDITIONAL_CERNLIBLIBS="-L${enableval}/lib -Wl,-rpath -Wl,${enableval}/lib -lpacklib_noshift -lmathlib -lkernlib_noshift -lXm"
                  cernlib=true;
                  AC_MSG_RESULT(${enableval});
		else
                CONDITIONAL_CERNLIBLIBS="${enableval}/lib/libpacklib_noshift.a ${enableval}/lib/libmathlib.a ${enableval}/lib/libkernlib_noshift.a"
                cernlib=true;
                AC_MSG_RESULT(${enableval});
		fi
              elif test -f "${enableval}/lib/libkernlib.a"; then
	        if test -f "${enableval}/lib/libkernlib.so"; then
                  CONDITIONAL_CERNLIBLIBS="-L${enableval}/lib -Wl,-rpath -Wl,${enableval}/lib -lpacklib -lmathlib -lkernlib -lXm"
                  cernlib=true;
                  AC_MSG_RESULT(${enableval});
		else
                CONDITIONAL_CERNLIBLIBS="${enableval}/lib/libpacklib.a ${enableval}/lib/libmathlib.a ${enableval}/lib/libkernlib.a"
                cernlib=true;
                AC_MSG_RESULT(${enableval});
		fi
              else
                AC_MSG_ERROR(Did not find '${enableval}/lib/libkernlib.a'.); 
              fi;
            else
              AC_MSG_ERROR(Did not find cernlib directory '${enableval}'.);
            fi;
      esac;
    ],
    [ cernlib=false ]
  )
  if test "$cernlib" = "true" ; then
    AC_DEFINE([USING__CERNLIB], "1", [cernlib found])
  fi
  AC_SUBST(CONDITIONAL_CERNLIBLIBS)
  AM_CONDITIONAL(CERNLIB_SUPPORT, test "$cernlib" = "true")

  AC_ARG_ENABLE(
    pgs,
    AS_HELP_STRING([--enable-pgs=/path/to/pgs],[Enable pgs.]),
    [ AC_MSG_CHECKING(for PGS installation directory);
      case "${enableval}" in
        no) AC_MSG_RESULT(PGS not enabled); pgs=false;;
        *)  if test -d "${enableval}"; then
              if test -f "${enableval}/lib/libstdhep.a"; then
                CONDITIONAL_PGSLIBS="${enableval}/lib/libstdhep.a ${enableval}/lib/libFmcfio.a"
                pgs=true;
                AC_MSG_RESULT(${enableval});
              else
                AC_MSG_ERROR(Did not find '${enableval}/lib/libstdhep.a'.); 
              fi;
            else
              AC_MSG_ERROR(Did not find PGS directory '${enableval}'.);
            fi;
      esac;
    ],
    [ pgs=false ]
  )
  AC_SUBST(CONDITIONAL_PGSLIBS)
  AM_CONDITIONAL(PGS_SUPPORT, test "$pgs" = "true")

  AC_ARG_ENABLE(
    delphes,
    AS_HELP_STRING([--enable-delphes=/path/to/delphes],[Enable delphes.]),
    [ AC_MSG_CHECKING(for DELPHES installation directory);
      case "${enableval}" in
        no) AC_MSG_RESULT(DELPHES not enabled); delphes=false;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_DELPHESLIBS="-Wl,-rpath -Wl,${enableval}/lib -L${enableval}/lib -lUtilities"
              CONDITIONAL_DELPHESINCS="-I${enableval}"
              delphes=true;
              AC_MSG_RESULT(${enableval});
            else
              AC_MSG_ERROR(Did not find DELPHES directory '${enableval}'.);
            fi;
      esac;
    ],
    [ delphes=false ]
  )
  if test "$delphes" = "true" ; then
    AC_DEFINE([USING__DELPHES], "1", [using delphes])
  fi
  AM_CONDITIONAL(DELPHES_SUPPORT, test "$delphes" = "true")
  AC_SUBST(CONDITIONAL_DELPHESLIBS)
  AC_SUBST(CONDITIONAL_DELPHESINCS)

  AC_ARG_ENABLE(
    gzip,
    AS_HELP_STRING([--enable-gzip],[Enable gzip support (for compressed event output)]),
    [ case "${enableval}" in
        no)   AC_MSG_RESULT(gzip not enabled); zlib=false ;;
        yes)  AC_CHECK_LIB(z, inflateEnd, [libz_found=yes], [libz_found=no])
              AC_CHECK_HEADER(zlib.h, [zlibh_found=yes], [zlibh_found=no])
              if test "$libz_found" = "yes" -a "$zlibh_found" = "yes"; then
                zlib=true;
                CONDITIONAL_GZIPLIBS="-lz";
              else
                AC_MSG_ERROR(Header zlib.h and/or library libz not found. Configure without --disable-gzip or install zlib (and its devel package, e.g. zlib-devel, zlib-dev or zlib1g-dev) if you want compressed output.);
              fi;;
	*) if test -d "${enableval}"; then
             ZLIB_OLD_LDFLAGS=$LDFLAGS; ZLIB_OLD_CPPFLAGS=$CPPFLAGS;
             LDFLAGS="$LDFLAGS -L${enableval}/lib"
             CPPFLAGS="$CPPFLAGS -I${enableval}/include"
             AC_CHECK_LIB(z,inflateEnd,zlib_cv_libz=yes,zlib_cv_libz=no)
             AC_CHECK_HEADER(zlib.h,zlib_cv_zlib_h=yes,zlib_cv_zlib_h=no)
             if test "$zlib_cv_libz" = "yes" && test "$zlib_cv_zlib_h" = "yes"
             then
               zlib=true;
               CONDITIONAL_GZIPLIBS="-Wl,-rpath -Wl,${enableval}/lib -L${enableval}/lib -lz"
               CONDITIONAL_GZIPINCS="-I${enableval}/include"
	       AC_MSG_RESULT(Using zlib from ${enableval})
             else
               AC_MSG_ERROR(Header zlib.h and/or library libz not found. Configure without --disable-gzip or install zlib.);
             fi
             LDFLAGS="$ZLIB_OLD_LDFLAGS"; CPPFLAGS="$ZLIB_OLD_CPPFLAGS";
           else
             AC_MSG_ERROR(no such directory '${enableval}');
           fi;;
      esac ],
    [ zlib=false ]
  )
  if test "$zlib" = "true" ; then
    AC_DEFINE([USING__GZIP], "1", [using gzip])
  fi
  AM_CONDITIONAL(GZIP_SUPPORT, test "$zlib" = "true")
  AC_SUBST(CONDITIONAL_GZIPLIBS)
  AC_SUBST(CONDITIONAL_GZIPINCS)

  AC_ARG_ENABLE(
    pythia,
    AS_HELP_STRING([--enable-pythia],[Enable fragmentation/decay interface to
    Pythia.]),
    [ AC_MSG_CHECKING(whether to enable Pythia interface);
      case "${enableval}" in
        no)   AC_MSG_RESULT(no); pythia=false ;;
        yes)  AC_MSG_RESULT(yes); pythia=true ;;
      esac ],
    [ pythia=false ]
  )
  if test "$pythia" = "true" ; then
    AC_DEFINE([USING__PYTHIA], "1", [Pythia interface enabled])
  fi
  AM_CONDITIONAL(PYTHIA_SUPPORT, test "$pythia" = "true")

  AC_ARG_ENABLE(
    pythia8,
    AS_HELP_STRING([--enable-pythia8=/path/to/pythia8], [Enable Pythia8 support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for Pythia8 installation directory);
      pythia8=true;
      case "${enableval}" in
        no)  AC_MSG_RESULT(Pythia8 not enabled); pythia8=false; pythia82=false; pythia83=false ;;
        yes) if test -x "`which pythia8-config`"; then
               CONDITIONAL_PYTHIA8DIR=`pythia8-config --prefix`;
             fi;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_PYTHIA8DIR=${enableval};
            fi;;
      esac;
      if test -n "$CONDITIONAL_PYTHIA8DIR"; then
        if test -x "$CONDITIONAL_PYTHIA8DIR/bin/pythia8-config"; then
          CONDITIONAL_PYTHIA8LDADD="$($CONDITIONAL_PYTHIA8DIR/bin/pythia8-config --ldflags)";
          CONDITIONAL_PYTHIA8CPPFLAGS="$($CONDITIONAL_PYTHIA8DIR/bin/pythia8-config --cxxflags)";
          AC_MSG_RESULT([${CONDITIONAL_PYTHIA8DIR}]);
        else
          AC_MSG_RESULT([Unable to find pythia8-config in the specifed path ${CONDITIONAL_PYTHIA8DIR}]);
          AC_MSG_CHECKING(Trying to proceed without pythia8-config);
          pythia8noconfiglibs=false;
          pythia8noconfigincludes=false;
          if test -f "${enableval}/lib/libpythia8.so"; then
             CONDITIONAL_PYTHIA8LDADD="-L${enableval}/lib -lpythia8";
             pythia8noconfiglibs=true;
          fi;
          if test -f "${enableval}/lib/libpythia8.dyld"; then
             CONDITIONAL_PYTHIA8LDADD="-L${enableval}/lib -lpythia8";
             pythia8noconfiglibs=true;
          fi;
          if test -f "${enableval}/lib64/libpythia8.so"; then
             CONDITIONAL_PYTHIA8LDADD="-L${enableval}/lib64 -lpythia8";
             pythia8noconfiglibs=true;
          fi;
          if test -f "${enableval}/include/Pythia8/Pythia.h"; then
             CONDITIONAL_PYTHIA8CPPFLAGS="-I${enableval}/Pythia8";
             pythia8noconfigincludes=true;
          fi;
          if "$pythia8noconfiglibs" = "true" &&  "$pythia8noconfigincludes" = "true"; then
             AC_MSG_RESULT([Found Pythia8 libraries and includes in ${CONDITIONAL_PYTHIA8DIR}]);
          else
             AC_MSG_ERROR(Unable to find Pythia8 headers and libraries from the specified path. );
          fi;
        fi;
      fi;
    ],
    [ pythia8=false ]
  )
  AC_SUBST(CONDITIONAL_PYTHIA8LDADD)
  AC_SUBST(CONDITIONAL_PYTHIA8CPPFLAGS)
  AM_CONDITIONAL(PYTHIA8_SUPPORT, test "$pythia8" = "true")
  if test "$pythia8" = "true" ; then
    AC_DEFINE([USING__PYTHIA8], "1", [Pythia8 interface enabled])
  fi

  AC_ARG_ENABLE(
    hepevtsize,
    AS_HELP_STRING([--enable-hepevtsize=HEPEVT_SIZE],[HEPEVT common block size @<:@default=10000@:>@]),
    [ AC_MSG_CHECKING(whether HEPEVT common block size is defined);
      if test ${enableval} -gt 0 2>/dev/null ; then
         HEPEVT_CB_SIZE=${enableval}
      	 AC_MSG_RESULT(${HEPEVT_CB_SIZE})
      fi
    ],
    [ HEPEVT_CB_SIZE=10000 ]
  )
  if test "x$HEPEVT_CB_SIZE" = "xno" ; then
        exit 1
  else
  	AC_DEFINE_UNQUOTED(HEPEVT_CB_SIZE, ${HEPEVT_CB_SIZE} , [HEPEVT common block size])
  fi
  AC_SUBST(HEPEVT_CB_SIZE)

  AC_ARG_ENABLE(
    binreloc,
    AS_HELP_STRING([--enable-binreloc],[Enable binrelocing]),
    [ AC_MSG_CHECKING(whether to install relocatable Sherpa)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); binreloc=false ;;
        yes) AC_MSG_RESULT(yes); binreloc=true ;;
      esac ],
    [ AC_MSG_CHECKING(whether to install relocatable Sherpa); AC_MSG_RESULT(no); binreloc=false ] 
  )
  if test "$binreloc" = "true" ; then
    AC_DEFINE([ENABLE_BINRELOC], "1", [binreloc activation])
  fi

  AC_ARG_ENABLE(pyext,
    AS_HELP_STRING([--enable-pyext],[Enable Python API]),
    [ AC_MSG_CHECKING(for Python extension)
      case "${enableval}" in
        no) AC_MSG_RESULT(no); pyext=false ;;
        yes)  AC_MSG_RESULT(yes); pyext=true ;;
      esac ],
    [ AC_MSG_CHECKING(for Python extension); AC_MSG_RESULT(no); pyext=false])
  if test x$pyext == xtrue; then
    AM_PATH_PYTHON
    AX_PYTHON_DEVEL
    dnl Note that we need swig 2.0.12 or later for C++11 compatibility
    AX_PKG_SWIG([2.0.12],[],[ AC_MSG_ERROR([SWIG is required to build..]) ])
    AX_SWIG_ENABLE_CXX
    AX_SWIG_MULTI_MODULE_SUPPORT
    AX_SWIG_PYTHON
  fi
  AM_CONDITIONAL(ENABLE_PYEXT, [test x$pyext == xtrue])

  AC_ARG_WITH([libzip],
    AS_HELP_STRING(
      [--with-libzip=@<:@ARG@:>@],
      [use libzip library from @<:@ARG@:>@]
    ),
    [
      if test "$withval" = "install"; then
      ac_libzip_path=$ac_default_prefix;
      test "x$prefix" != xNONE && ac_libzip_path=$prefix;
      if ! test -f ${ac_libzip_path}/include/zip.h; then
        wget --no-check-certificate https://libzip.org/download/libzip-1.2.0.tar.gz
        tar xzf libzip-1.2.0.tar.gz;
        cd libzip-1.2.0;
        ./configure --prefix=${ac_libzip_path} || exit;
        make || exit; make install || exit;
        mv ${ac_libzip_path}/lib/libzip/include/zipconf.h ${ac_libzip_path}/include/;
        cd ..;
        rm -rf libzip-1.2.0.tar.gz libzip-1.2.0;
        echo "Successfully installed libzip into ${ac_libzip_path}."
      fi
      else
        ac_libzip_path="$withval"
      fi
    ],
    [ ac_libzip_path=/usr; ]
  )
  if ! test -f ${ac_libzip_path}/include/zip.h; then
    AC_MSG_ERROR(Did not find required dependency libzip in ${ac_libzip_path}. Please specify its installation prefix using '--with-libzip=/path' or enable its automatic installation using '--with-libzip=install'.)
  fi
  LIBZIP_CPPFLAGS="-I$ac_libzip_path/include -I$ac_libzip_path/lib/libzip/include"
  LIBZIP_LDFLAGS="-L$ac_libzip_path/lib -L$ac_libzip_path/lib64 -lzip"
  AC_SUBST(LIBZIP_CPPFLAGS)
  AC_SUBST(LIBZIP_LDFLAGS)

])
