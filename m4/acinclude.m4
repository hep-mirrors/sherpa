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
      AC_DEFINE([ARCH_DARWIN], "1", [Architecture identified as Darwin MacOS])
      AC_DEFINE([LIB_SUFFIX], ".dylib", [library suffix set to .dylib]) 
      AC_DEFINE([LD_PATH_NAME], "DYLD_LIBRARY_PATH", [ld path name set to DYLD_LIBRARY_PATH]) ;;
    *linux*:*:*)
      echo "checking for architecture...  Linux"
      if test "x$LDFLAGS" = "x"; then
        AM_LDFLAGS="-rdynamic -Wl,--no-as-needed"
      fi
      SEDCOMMAND="sed -i -r"
      AC_DEFINE([ARCH_LINUX], "1", [Architecture identified as Linux])
      AC_DEFINE([LIB_SUFFIX], ".so", [library suffix set to .so]) 
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
    *)
      echo "checking for architecture...  unknown"
      echo "hosts system type $build not yet supported, assuming unix behaviour."
      echo "possible failure due to unknown compiler/linker characteristics."
      echo "please inform us about build results at sherpa@projects.hepforge.org"
      echo "(will continue in 10 seconds)"
      sleep 10
      if test "x$LDFLAGS" = "x"; then
        AM_LDFLAGS="-rdynamic -Wl,--no-as-needed"
      fi
      SEDCOMMAND="sed -i -r"
      AC_DEFINE([ARCH_UNIX], "1", [Architecture identified as Unix])
      AC_DEFINE([LIB_SUFFIX], ".so", [library suffix set to .so]) 
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
  esac
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
	\${AMEGICBUILDDIR}/Amplitude/Zfunctions/libZfunctions.la \
	\${AMEGICBUILDDIR}/Cluster/libAmegicCluster.la"
  AC_SUBST(AMEGICDIR)
  AC_SUBST(AMEGICBUILDDIR)
  AC_SUBST(AMEGICLIBS)

  AMISICDIR="\${top_srcdir}/AMISIC++"
  AMISICBUILDDIR="\${top_builddir}/AMISIC++"
  AMISICLIBS="\${AMISICBUILDDIR}/Main/libAmisic.la \
	\${AMISICBUILDDIR}/Tools/libAmisicTools.la \
	\${AMISICBUILDDIR}/Model/libAmisicModel.la"
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
	\${ATOOLSBUILDDIR}/Org/libToolsOrg.la"
  AC_SUBST(ATOOLSDIR)
  AC_SUBST(ATOOLSBUILDDIR)
  AC_SUBST(ATOOLSLIBS)
  
  BEAMDIR="\${top_srcdir}/BEAM"
  BEAMBUILDDIR="\${top_builddir}/BEAM"
  BEAMLIBS="\${BEAMBUILDDIR}/Main/libBeam.la"
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
	\${EXTRAXSBUILDDIR}/Cluster/libExtraXSCluster.la \
	\${EXTRAXSBUILDDIR}/NLO/libExtraXSNLO.la"
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
  

  COMIXDIR="\${top_srcdir}/COMIX"
  COMIXBUILDDIR="\${top_builddir}/COMIX"
  COMIXLIBS="\${COMIXBUILDDIR}/Amplitude/libComixAmplitude.la \
	\${COMIXBUILDDIR}/Phasespace/libComixPhasespace.la \
	\${COMIXBUILDDIR}/Main/libComix.la \
	\${COMIXBUILDDIR}/Cluster/libComixCluster.la"
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
  MODELLIBS="\${MODELBUILDDIR}/Main/libModelMain.la \
	\${MODELBUILDDIR}/UFO/libModelUFO.la"
  AC_SUBST(MODELDIR)
  AC_SUBST(MODELBUILDDIR)
  AC_SUBST(MODELLIBS)
  
  PDFDIR="\${top_srcdir}/PDF"
  PDFBUILDDIR="\${top_builddir}/PDF"
  PDFINCS="-I\${PDFDIR}/Main -I\${PDFDIR}/Remnant"
  PDFLIBS="\${PDFBUILDDIR}/Main/libPDF.la \
	\${PDFBUILDDIR}/Remnant/libRemnant.la"
  AC_SUBST(PDFDIR)
  AC_SUBST(PDFBUILDDIR)
  AC_SUBST(PDFLIBS)
  
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

  AC_DEFINE_UNQUOTED([SHERPA_VERSION], ["`echo AC_PACKAGE_VERSION | cut -d. -f1`"], [Sherpa version])
  AC_DEFINE_UNQUOTED([SHERPA_SUBVERSION], ["`echo AC_PACKAGE_VERSION | cut -d. -f2,3`"], [Sherpa subversion])
  AC_DEFINE_UNQUOTED([SHERPA_PREFIX], "$SHERPAPREFIX", [Sherpa installation prefix])
  AC_DEFINE_UNQUOTED([SHERPA_INCLUDE_PATH], "$INCLUDEDIR", [Sherpa include directory])
  AC_DEFINE_UNQUOTED([SHERPA_LIBRARY_PATH], "$LIBDIR", [Sherpa library directory])
  AC_DEFINE_UNQUOTED([SHERPA_SHARE_PATH], "$DATADIR", [Sherpa data directory])
  AC_DEFINE([USING__COLOUR], "1", [Using colour])
  
  AM_CPPFLAGS="-I\$(top_srcdir)"
  AC_SUBST(AM_CPPFLAGS)

  AM_CXXFLAGS="-g -O2"
  AC_SUBST(AM_CXXFLAGS)

  localincdir="\$(pkgincludedir)/\$(subdir)"
  AC_SUBST(localincdir)
])



dnl Conditional compiling and linking

AC_DEFUN([SHERPA_SETUP_CONFIGURE_OPTIONS],
[
  AC_ARG_ENABLE(
    versioning,
    AC_HELP_STRING([--enable-versioning], [Add version tag to executables and library/header directories, such that multiple Sherpa versions can live in the same prefix.]),
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
  AC_SUBST(VERSIONING)

  AC_ARG_ENABLE(
    multithread,
    AC_HELP_STRING([--enable-multithread], [Enable multithreading]),
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
    AC_HELP_STRING([--enable-analysis], [Enable analysis]),
    [ AC_MSG_CHECKING(for analysis)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); analysis=false ;;
        yes) AC_MSG_RESULT(yes); analysis=true ;;
      esac ],
    [ AC_MSG_CHECKING(for analysis); AC_MSG_RESULT(no); analysis=false ]
  )
  AM_CONDITIONAL(USING__Analysis, test "$analysis" = "true" )

  AC_ARG_ENABLE(
    hepmc2,
    AC_HELP_STRING([--enable-hepmc2=/path/to/hepmc], [Enable HepMC (version 2.x) support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for HepMC2 installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(HepMC2 not enabled); hepmc2=false ;;
        yes)  if test -d "$HEPMC2DIR"; then
                CONDITIONAL_HEPMC2DIR="$HEPMC2DIR"
                CONDITIONAL_HEPMC2INCS="-I$HEPMC2DIR/include"
                CONDITIONAL_HEPMC2LIBS="-L$HEPMC2DIR/lib -R$HEPMC2DIR/lib -L$HEPMC2DIR/lib64 -R$HEPMC2DIR/lib64 -lHepMC";
              else
                AC_MSG_ERROR(\$HEPMC2DIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_HEPMC2DIR}]); hepmc2=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_HEPMC2DIR="${enableval}"
                CONDITIONAL_HEPMC2INCS="-I${enableval}/include"
                CONDITIONAL_HEPMC2LIBS="-L${enableval}/lib -R${enableval}/lib -L${enableval}/lib64 -R${enableval}/lib64 -lHepMC";
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
    rivet,
    AC_HELP_STRING([--enable-rivet=/path/to/rivet], [Enable Rivet support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for Rivet installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(Rivet not enabled); rivet=false ;;
        yes) if test -x "`which rivet-config`"; then
               CONDITIONAL_RIVETDIR=`rivet-config --prefix`;
             fi;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_RIVETDIR=${enableval};
            fi;;
      esac;
      if test -x "$CONDITIONAL_RIVETDIR/bin/rivet-config"; then
        CONDITIONAL_RIVETLDADD="$($CONDITIONAL_RIVETDIR/bin/rivet-config --ldflags) $($CONDITIONAL_RIVETDIR/bin/rivet-config --ldadd)";
        CONDITIONAL_RIVETCPPFLAGS="$($CONDITIONAL_RIVETDIR/bin/rivet-config --cppflags)";
        AC_MSG_RESULT([${CONDITIONAL_RIVETDIR}]); rivet=true;
	"$CONDITIONAL_RIVETDIR/bin/rivet-config" --version | grep -q '^1\.' || rivetyoda=true
      else
        AC_MSG_ERROR(Unable to use Rivet from specified path.);
      fi;
      rivetincludedir=$($CONDITIONAL_RIVETDIR/bin/rivet-config --includedir)
      if grep -q -s setIgnoreBeams $rivetincludedir/Rivet/AnalysisHandler.hh; then
        rivetsetsow=true;
      fi
    ],
    [ rivet=false ]
  )
  if test "$rivet" = "true" ; then
    AC_DEFINE([USING__RIVET], "1", [using Rivet])
  fi
  if test "$rivetsetsow" = "true" ; then
    AC_DEFINE([USING__RIVET__SETSOW], "1", [setSumOfWeights function available in Rivet])
  fi
  if test "$rivetyoda" = "true" ; then
    AC_DEFINE([USING__RIVET__YODA], "1", [Rivet uses YODA as its histogramming backend])
  fi
  AC_SUBST(CONDITIONAL_RIVETLDADD)
  AC_SUBST(CONDITIONAL_RIVETCPPFLAGS)
  AM_CONDITIONAL(RIVET_SUPPORT, test "$rivet" = "true")
  

  AC_ARG_ENABLE(
    fastjet,
    AC_HELP_STRING([--enable-fastjet=/path/to/fastjet], [Enable FASTJET.]),
    [ AC_MSG_CHECKING(for FASTJET installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(FASTJET not enabled); fastjet=false ;;
        yes)  if test -d "$FASTJETDIR"; then
                CONDITIONAL_FASTJETDIR="$FASTJETDIR"
                CONDITIONAL_FASTJETINCS="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --cxxflags)";
                CONDITIONAL_FASTJETLIBS="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --libs --plugins=yes)"
                CONDITIONAL_FASTJETVERSION="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --version)";
              else
                AC_MSG_ERROR(\$FASTJETDIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_FASTJETDIR}]); fastjet=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_FASTJETDIR="${enableval}"
                CONDITIONAL_FASTJETINCS="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --cxxflags)";
                CONDITIONAL_FASTJETLIBS="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --libs --plugins=yes)"
                CONDITIONAL_FASTJETVERSION="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --version)";
              else
                AC_MSG_ERROR(${enableval} is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_FASTJETDIR}]); fastjet=true;;
      esac
      ],
    [ fastjet=false ]
  )
  if test "$fastjet" = "true" ; then
    AC_DEFINE([USING__FASTJET], "1", [Using FASTJET])
  fi
  if test "$(echo $CONDITIONAL_FASTJETVERSION | cut -d . -f 1)" = "3"; then
    AC_DEFINE([USING__FASTJET__3], "1", [Using FASTJET 3])
  fi
  AC_SUBST(CONDITIONAL_FASTJETDIR)
  AC_SUBST(CONDITIONAL_FASTJETINCS)
  AC_SUBST(CONDITIONAL_FASTJETLIBS)
  AM_CONDITIONAL(FASTJET_SUPPORT, test "$fastjet" = "true")

  AC_ARG_ENABLE(
    blackhat,
    AC_HELP_STRING([--enable-blackhat=/path/to/blackhat], [Enable BLACKHAT.]),
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
    AC_HELP_STRING([--enable-openloops=/path/to/openloops], [Enable OpenLoops.]),
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
    mcfm,
    AC_HELP_STRING([--enable-mcfm=/path/to/mcfm], [Enable MCFM.]),
    [ AC_MSG_CHECKING(for MCFM installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(MCFM not enabled); mcfm=false ;;
        yes)  if test -d "$MCFMDIR"; then
                CONDITIONAL_MCFMDIR="$MCFMDIR"
                CONDITIONAL_MCFMLIBS="$CONDITIONAL_MCFMDIR/lib/libMCFM.a"
              else
                AC_MSG_ERROR(\$MCFMDIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_MCFMDIR}]); mcfm=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_MCFMDIR="${enableval}"
                CONDITIONAL_MCFMLIBS="$CONDITIONAL_MCFMDIR/lib/libMCFM.a"
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
  fi
  AC_SUBST(CONDITIONAL_MCFMDIR)
  AC_SUBST(CONDITIONAL_MCFMLIBS)
  AM_CONDITIONAL(MCFM_SUPPORT, test "$mcfm" = "true")

  AC_ARG_ENABLE(
    lhole,
    AC_HELP_STRING([--enable-lhole], [Enable Les Houches One-Loop Generator interface.]),
    [ AC_MSG_CHECKING(for LHOLE)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); lhole=false ;;
        yes) AC_MSG_RESULT(yes); lhole=true ;;
      esac ],
    [ AC_MSG_CHECKING(for LHOLE); AC_MSG_RESULT(no); lhole=false ]
  )
  AM_CONDITIONAL(USING__LHOLE, test "$lhole" = "true" )

  AC_ARG_ENABLE(
    root,
    AC_HELP_STRING([--enable-root=/path/to/root], [Enable ROOT support and specify where it is installed if non-standard.]),
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
  

  lhapdfversion=5
  AC_ARG_ENABLE(
    lhapdf,
    AC_HELP_STRING([--enable-lhapdf=/path/to/lhapdf], [Enable LHAPDF support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for LHAPDF installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(LHAPDF not enabled); lhapdf=false ;;
        yes) if test -x "`which lhapdf-config`"; then
               CONDITIONAL_LHAPDFDIR=`lhapdf-config --prefix`;
             fi;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_LHAPDFDIR=${enableval};
            fi;;
      esac;

      if test -x "$CONDITIONAL_LHAPDFDIR/bin/lhapdf-config"; then
        CONDITIONAL_LHAPDFLIBS="$($CONDITIONAL_LHAPDFDIR/bin/lhapdf-config --ldflags)";
        CONDITIONAL_LHAPDFINCS="$($CONDITIONAL_LHAPDFDIR/bin/lhapdf-config --cppflags)";
        lhapdfversion="$($CONDITIONAL_LHAPDFDIR/bin/lhapdf-config --version)";
        lhapdfversion=${lhapdfversion:0:1}
        AC_MSG_RESULT([${CONDITIONAL_LHAPDFDIR}]); lhapdf=true;
      else
        AC_MSG_ERROR(Unable to use LHAPDF from specified path.);
      fi;
    ],
    [ lhapdf=false ]
  )
  if test "$lhapdf" = "true" ; then
    AC_DEFINE_UNQUOTED([LHAPDF_PATH], "$CONDITIONAL_LHAPDFDIR", [LHAPDF directory])
    AC_DEFINE([USING__LHAPDF], "1", [using LHAPDF])
    if test [ "$lhapdfversion" -ge "6" ] ; then
      AC_DEFINE(USING__LHAPDF6, "1", [using LHAPDF6])
    fi
  fi
  AC_SUBST(CONDITIONAL_LHAPDFDIR)
  AC_SUBST(CONDITIONAL_LHAPDFLIBS)
  AC_SUBST(CONDITIONAL_LHAPDFINCS)
  AM_CONDITIONAL(LHAPDF_SUPPORT, test "$lhapdf" = "true")
  AM_CONDITIONAL(LHAPDF6_SUPPORT, test [ "$lhapdfversion" -ge "6" ])

  AC_ARG_ENABLE(
    hztool,
    AC_HELP_STRING([--enable-hztool=/path/to/hztool], [Enable hztool for analysis.]),
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
                AC_MSG_ERROR(Did not find '${enableval}/libhztool.so'.); 
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
    AC_HELP_STRING([--enable-cernlib=/path/to/cernlib], [Enable cernlib.]),
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
    AC_HELP_STRING([--enable-pgs=/path/to/pgs], [Enable pgs.]),
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
    AC_HELP_STRING([--enable-delphes=/path/to/delphes], [Enable delphes.]),
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
    AC_HELP_STRING([--enable-gzip], [Enable gzip support (for compressed event output)]),
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
    AC_HELP_STRING([--enable-pythia], [Enable fragmentation/decay interface to
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
    hepevtsize,
    AC_HELP_STRING([--enable-hepevtsize=HEPEVT_SIZE], [HEPEVT common block size @<:@default=10000@:>@]),
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
    AC_HELP_STRING([--enable-binreloc], [Enable binrelocing]),
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

])
