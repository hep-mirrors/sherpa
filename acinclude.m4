dnl workaround for old automake on darwin

AC_DEFUN([AM_CONFIG_HEADERS], [AC_CONFIG_HEADERS($@)])

dnl set flags according to build environment

AC_DEFUN([SHERPA_SETUP_BUILDSYSTEM],
[
  case "$build_os:$build_cpu:$build_vendor" in
    *darwin*:*:*)
      echo "checking for architecture... Darwin MacOS"
      AM_LDFLAGS=""
      ldflags="-dynamic -flat_namespace"
      SEDCOMMAND="sed -i.bak -E"
      AC_DEFINE([ARCH_DARWIN], "1", [Architecture identified as Darwin MacOS])
      AC_DEFINE([LIB_SUFFIX], ".dylib", [library suffix set to .dylib]) 
      AC_DEFINE([LD_PATH_NAME], "DYLD_LIBRARY_PATH", [ld path name set to DYLD_LIBRARY_PATH]) ;;
    *linux*:*:*)
      echo "checking for architecture...  Linux"
      AM_LDFLAGS="-rdynamic -Xlinker --no-as-needed -Xlinker --add-needed"
      ldflags="-rdynamic"
      SEDCOMMAND="sed -i -r"
      AC_DEFINE([ARCH_LINUX], "1", [Architecture identified as Linux])
      AC_DEFINE([LIB_SUFFIX], ".so", [library suffix set to .so]) 
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
    *)
      echo "checking for architecture...  unknown"
      echo "hosts system type $build not yet supported, assuming unix behaviour."
      echo "possible failure due to unknown compiler/linker characteristics."
      echo "please inform us about build results at info@sherpa-mc.de"
      echo "(will continue in 10 seconds)"
      sleep 10
      AM_LDFLAGS="-rdynamic -Xlinker --no-as-needed -Xlinker --add-needed"
      ldflags="-rdynamic"
      SEDCOMMAND="sed -i -r"
      AC_DEFINE([ARCH_UNIX], "1", [Architecture identified as Unix])
      AC_DEFINE([LIB_SUFFIX], ".so", [library suffix set to .so]) 
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
  esac
  AC_SUBST(AM_LDFLAGS)
  AC_SUBST(ldflags)
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
  AMEGICLIBS="-L\${AMEGICBUILDDIR}/Main -L\${AMEGICBUILDDIR}/DipoleSubtraction \ 
	      -L\${AMEGICBUILDDIR}/Amplitude -L\${AMEGICBUILDDIR}/Phasespace \
              -L\${AMEGICBUILDDIR}/String -L\${AMEGICBUILDDIR}/Amplitude/Zfunctions -L\${AMEGICBUILDDIR}/Cluster \
              -lAmegic -lDipoleSubtraction -lAmplitude -lAmegicPSGen -lZfunctions -lString -lAmegicCluster"
  AC_SUBST(AMEGICDIR)
  AC_SUBST(AMEGICBUILDDIR)
  AC_SUBST(AMEGICLIBS)

  AMISICDIR="\${top_srcdir}/AMISIC++"
  AMISICBUILDDIR="\${top_builddir}/AMISIC++"
  AMISICLIBS="-L\${AMISICBUILDDIR}/Main -L\${AMISICBUILDDIR}/Tools -L\${AMISICBUILDDIR}/Model \
              -lAmisic -lAmisicModel -lAmisicTools"
  AC_SUBST(AMISICDIR)
  AC_SUBST(AMISICBUILDDIR)
  AC_SUBST(AMISICLIBS)

  AHADICDIR="\${top_srcdir}/AHADIC++"
  AHADICBUILDDIR="\${top_builddir}/AHADIC++"
  AHADICLIBS="-L\${AHADICBUILDDIR}/Main -L\${AHADICBUILDDIR}/Tools -L\${AHADICBUILDDIR}/Formation -L\${AHADICBUILDDIR}/Decays \
              -lAhadicMain -lAhadicTools -lAhadicFormation -lAhadicDecays"
  AC_SUBST(AHADICDIR)
  AC_SUBST(AHADICBUILDDIR)
  AC_SUBST(AHADICLIBS)
  
  ATOOLSDIR="\${top_srcdir}/ATOOLS"
  ATOOLSBUILDDIR="\${top_builddir}/ATOOLS"
  ATOOLSLIBS="-L\${ATOOLSBUILDDIR}/Phys -L\${ATOOLSBUILDDIR}/Math \
              -L\${ATOOLSBUILDDIR}/Org \
              -lToolsPhys -lToolsMath -lToolsOrg"
  AC_SUBST(ATOOLSDIR)
  AC_SUBST(ATOOLSBUILDDIR)
  AC_SUBST(ATOOLSLIBS)
  
  BEAMDIR="\${top_srcdir}/BEAM"
  BEAMBUILDDIR="\${top_builddir}/BEAM"
  BEAMLIBS="-L\${BEAMBUILDDIR}/Main -lBeam"
  AC_SUBST(BEAMDIR)
  AC_SUBST(BEAMBUILDDIR)
  AC_SUBST(BEAMLIBS)

  METOOLSDIR="\${top_srcdir}/METOOLS"
  METOOLSBUILDDIR="\${top_builddir}/METOOLS"
  METOOLSLIBS="-L\${METOOLSBUILDDIR}/Explicit -lMEToolsExplicit \
                  -L\${METOOLSBUILDDIR}/Currents -lMEToolsCurrents \
                  -L\${METOOLSBUILDDIR}/Vertices -lMEToolsVertices \
                  -L\${METOOLSBUILDDIR}/Colors  -lMEToolsColors \
                  -L\${METOOLSBUILDDIR}/SpinCorrelations -lMEToolsSpinCorrelations \
                  -L\${METOOLSBUILDDIR}/Loops -lMEToolsLoops \
                  -L\${METOOLSBUILDDIR}/Main -lMEToolsMain"
  AC_SUBST(METOOLSDIR)
  AC_SUBST(METOOLSBUILDDIR)
  AC_SUBST(METOOLSLIBS)
  
  EXTRAXSDIR="\${top_srcdir}/EXTRA_XS"
  EXTRAXSBUILDDIR="\${top_builddir}/EXTRA_XS"
  EXTRAXSLIBS="-L\${EXTRAXSBUILDDIR}/Main -lExtraXS \
               -L\${EXTRAXSBUILDDIR}/Two2Two -lExtraXS2_2 \
               -L\${EXTRAXSBUILDDIR}/One2Two -lExtraXS1_2 \
               -L\${EXTRAXSBUILDDIR}/One2Three -lExtraXS1_3 \
               -L\${EXTRAXSBUILDDIR}/Cluster -lExtraXSCluster \
               -L\${EXTRAXSBUILDDIR}/NLO -lExtraXSNLO"
  AC_SUBST(EXTRAXSDIR)
  AC_SUBST(EXTRAXSBUILDDIR)
  AC_SUBST(EXTRAXSLIBS)
  
  POWHEGDIR="\${top_srcdir}/POWHEG"
  POWHEGBUILDDIR="\${top_builddir}/POWHEG"
  POWHEGLIBS="-L\${POWHEGBUILDDIR}/Main -L\${POWHEGBUILDDIR}/Calculators -L\${POWHEGBUILDDIR}/Showers -L\${POWHEGBUILDDIR}/Tools \
		-lPOWHEGTools -lPOWHEGCalculators -lPOWHEGShowers -lPOWHEGMain"
  AC_SUBST(POWHEGDIR)
  AC_SUBST(POWHEGBUILDDIR)
  AC_SUBST(POWHEGLIBS)
  
  CSSDIR="\${top_srcdir}/CSSHOWER++"
  CSSBUILDDIR="\${top_builddir}/CSSHOWER++"
  CSSLIBS="-L\${CSSBUILDDIR}/Main -L\${CSSBUILDDIR}/Calculators -L\${CSSBUILDDIR}/Showers -L\${CSSBUILDDIR}/Tools \
		-lCSTools -lCSCalculators -lCSShowers -lCSMain"
  AC_SUBST(CSSDIR)
  AC_SUBST(CSSBUILDDIR)
  AC_SUBST(CSSLIBS)
  

  COMIXDIR="\${top_srcdir}/COMIX"
  COMIXBUILDDIR="\${top_builddir}/COMIX"
  COMIXLIBS="-L\${COMIXBUILDDIR}/Amplitude -L\${COMIXBUILDDIR}/Phasespace -L\${COMIXBUILDDIR}/Main -L\${COMIXBUILDDIR}/Cluster -lComixAmplitude -lComixPhasespace -lComix -lComixCluster"
  AC_SUBST(COMIXDIR)
  AC_SUBST(COMIXBUILDDIR)
  AC_SUBST(COMIXLIBS)
  
  HADRONSDIR="\${top_srcdir}/HADRONS++"
  HADRONSBUILDDIR="\${top_builddir}/HADRONS++"
  HADRONSLIBS="-L\${HADRONSBUILDDIR}/Main -L\${HADRONSBUILDDIR}/ME_Library \
               -L\${HADRONSBUILDDIR}/Current_Library -L\${HADRONSBUILDDIR}/PS_Library \
               -lHadronsMain -lHadronsMEs -lHadronsCurrents -lHadronsPSs"
  AC_SUBST(HADRONSDIR)
  AC_SUBST(HADRONSBUILDDIR)
  AC_SUBST(HADRONSLIBS)
  
  PHOTONSDIR="\${top_srcdir}/PHOTONS++"
  PHOTONSBUILDDIR="\${top_builddir}/PHOTONS++"
  PHOTONSLIBS="-L\${PHOTONSBUILDDIR}/Main -L\${PHOTONSBUILDDIR}/Tools \
               -L\${PHOTONSBUILDDIR}/PhaseSpace -L\${PHOTONSBUILDDIR}/MEs \
               -lPhotonsMain -lPhotonsTools -lPhotonsPhaseSpace -lPhotonsMEs"
  AC_SUBST(PHOTONSDIR)
  AC_SUBST(PHOTONSBUILDDIR)
  AC_SUBST(PHOTONSLIBS)
  
  MODELDIR="\${top_srcdir}/MODEL"
  MODELBUILDDIR="\${top_builddir}/MODEL"
  MODELLIBS="-L\${MODELBUILDDIR}/Main -L\${MODELBUILDDIR}/Interaction_Models \	
             -lModelMain -lModelInteractions"
  AC_SUBST(MODELDIR)
  AC_SUBST(MODELBUILDDIR)
  AC_SUBST(MODELLIBS)
  
  PDFDIR="\${top_srcdir}/PDF"
  PDFBUILDDIR="\${top_builddir}/PDF"
  PDFINCS="-I\${PDFDIR}/Main -I\${PDFDIR}/Remnant"
  PDFLIBS="-L\${PDFBUILDDIR}/Main -L\${PDFBUILDDIR}/Remnant \
           -lPDF -lRemnant"
  AC_SUBST(PDFDIR)
  AC_SUBST(PDFBUILDDIR)
  AC_SUBST(PDFLIBS)
  
  PHASICDIR="\${top_srcdir}/PHASIC++"
  PHASICBUILDDIR="\${top_builddir}/PHASIC++"
  PHASICLIBS="-L\${PHASICBUILDDIR}/Main -L\${PHASICBUILDDIR}/Channels \
	-L\${PHASICBUILDDIR}/Process -L\${PHASICBUILDDIR}/Scales \
	-L\${PHASICBUILDDIR}/Selectors \
	-lPhasicChannels -lPhasicMain -lPhasicProcess -lPhasicScales \
	-lPhasicSelectors \
    -L\${PHASICBUILDDIR}/Decays -lPhasicDecays"
  AC_SUBST(PHASICDIR)
  AC_SUBST(PHASICBUILDDIR)
  AC_SUBST(PHASICLIBS)
  
  SHERPADIR="\${top_srcdir}/SHERPA"
  SHERPABUILDDIR="\${top_builddir}/SHERPA"
  SHERPALIBS="-L\${SHERPABUILDDIR}/Single_Events -L\${SHERPABUILDDIR}/PerturbativePhysics \
              -L\${SHERPABUILDDIR}/LundTools -L\${SHERPABUILDDIR}/Tools -L\${SHERPABUILDDIR}/Main \
              -L\${SHERPABUILDDIR}/Initialization -L\${SHERPABUILDDIR}/SoftPhysics -L\${SHERPABUILDDIR}/HerwigTools \
              -lSherpaMain -lSherpaInitialization -lSherpaSingleEvents \
              -lSherpaPerturbativePhysics -lSherpaSoftPhysics -lLundTools -lSherpaTools"
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

  EXTRA_DIST="\$(SVNINCLUDE)"
  AC_SUBST(EXTRA_DIST)

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
    svninclude,
    AC_HELP_STRING([--enable-svninclude], [Add .svn directories to dist tarball.]),
    [ AC_MSG_CHECKING(whether to enable SVN synchronization)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no);
              SVNINCLUDE="";;
        yes) AC_MSG_RESULT(yes);
              SVNINCLUDE=".svn";;
      esac ],
    [ AC_MSG_CHECKING(whether to enable SVN synchronization); AC_MSG_RESULT(no); SVNINCLUDE="" ] 
  )
  AC_SUBST(SVNINCLUDE)

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
                CONDITIONAL_HEPMC2LIBS="-L$HEPMC2DIR/lib -R$HEPMC2DIR/lib -lHepMC";
              else
                AC_MSG_ERROR(\$HEPMC2DIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_HEPMC2DIR}]); hepmc2=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_HEPMC2DIR="${enableval}"
                CONDITIONAL_HEPMC2INCS="-I${enableval}/include"
                CONDITIONAL_HEPMC2LIBS="-L${enableval}/lib -R${enableval}/lib -lHepMC";
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
              else
                AC_MSG_ERROR(\$FASTJETDIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_FASTJETDIR}]); fastjet=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_FASTJETDIR="${enableval}"
                CONDITIONAL_FASTJETINCS="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --cxxflags)";
                CONDITIONAL_FASTJETLIBS="$($CONDITIONAL_FASTJETDIR/bin/fastjet-config --libs --plugins=yes)"
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
    phox,
    AC_HELP_STRING([--enable-phox=/path/to/phox], [Enable Phox.]),
    [ AC_MSG_CHECKING(for Phox installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(Phox not enabled); phox=false ;;
        yes)  if test -d "$PhoxDIR"; then
                CONDITIONAL_PhoxDIR="$PhoxDIR"
                CONDITIONAL_PhoxLIBS="$CONDITIONAL_PhoxDIR/lib/libPhox.a"
              else
                AC_MSG_ERROR(\$PhoxDIR is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_PhoxDIR}]); phox=true;;
        *)    if test -d "${enableval}"; then
                CONDITIONAL_PhoxDIR="${enableval}"
              else
                AC_MSG_ERROR(${enableval} is not a valid path.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_PhoxDIR}]); phox=true;;
      esac
      ],
    [ phox=false ]
  )
  if test "$phox" = "true" ; then
    AC_DEFINE([USING__Phox], "1", [Using Phox])
  fi
  AC_SUBST(CONDITIONAL_PhoxDIR)
  AC_SUBST(CONDITIONAL_PhoxLIBS)
  AM_CONDITIONAL(Phox_SUPPORT, test "$phox" = "true")

  AC_ARG_ENABLE(
    root,
    AC_HELP_STRING([--enable-root\[=/path/to/root\]], [Enable ROOT support and specify where it is installed if non-standard.]),
    [ AC_MSG_CHECKING(for ROOT installation directory)
      case "${enableval}" in
        no)  AC_MSG_RESULT(ROOT not enabled); root=false;;
        yes) if test -d "$ROOTSYS"; then
               CONDITIONAL_ROOTDIR=$ROOTSYS
               CONDITIONAL_ROOTINCS="-I$ROOTSYS/include -I$($ROOTSYS/bin/root-config --incdir)";
               CONDITIONAL_ROOTLIBS="-L$ROOTSYS/lib $($ROOTSYS/bin/root-config --glibs)"
               CONDITIONAL_ROOTFLAGS=-Wno-long-long
             elif test -x "`which root-config`"; then
               CONDITIONAL_ROOTDIR=`root-config --prefix`;
               CONDITIONAL_ROOTINCS=-I`root-config --incdir`;
               CONDITIONAL_ROOTLIBS=`root-config --glibs`;
               CONDITIONAL_ROOTFLAGS=-Wno-long-long
                if ! test -d "$CONDITIONAL_ROOTDIR"; then
                  AC_MSG_ERROR(root-config --prefix returned a path that is not available. Please check your ROOT installation and set \$ROOTSYS manually.);
                fi
             else
               AC_MSG_ERROR(\$ROOTSYS is not a valid path and root-config was not found.);
             fi;
             AC_MSG_RESULT([${CONDITIONAL_ROOTDIR}]); root=true;;
        *)   if test -d "${enableval}"; then
               CONDITIONAL_ROOTDIR="${enableval}"
               CONDITIONAL_ROOTINCS="-I${enableval}/include";
               CONDITIONAL_ROOTLIBS="-L${enableval}/lib $(${enableval}/bin/root-config --glibs)";
               CONDITIONAL_ROOTFLAGS="-Wno-long-long"
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
    lhapdf,
    AC_HELP_STRING([--enable-lhapdf=/path/to/lhapdf], [Enable LHAPDF support and specify where it is installed.]),
    [ AC_MSG_CHECKING(for LHAPDF installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(LHAPDF not enabled); lhapdf=false ;;
        yes) if test -d "$LHAPDFDIR"; then
               CONDITIONAL_LHAPDFDIR=$LHAPDFDIR;
             elif test -x "`which lhapdf-config`"; then
               CONDITIONAL_LHAPDFDIR=`lhapdf-config --prefix`;
             fi
             AC_MSG_RESULT([${CONDITIONAL_LHAPDFDIR}]); lhapdf=true;;
        *)  if test -d "${enableval}"; then
              CONDITIONAL_LHAPDFDIR=${enableval};
            fi;
      esac;
      if test -d "$CONDITIONAL_LHAPDFDIR"; then
        if test -f "$CONDITIONAL_LHAPDFDIR/lib/libLHAPDF.la"; then
          CONDITIONAL_LHAPDFLIBS="-L\${PDFBUILDDIR}/LHAPDF -lLHAPDFSherpa -L$CONDITIONAL_LHAPDFDIR/lib -lLHAPDF";
          CONDITIONAL_LHAPDFINCS="-I$CONDITIONAL_LHAPDFDIR/include";
        else
          CONDITIONAL_LHAPDFLIBS="-L\${PDFBUILDDIR}/LHAPDF -lLHAPDFSherpa $CONDITIONAL_LHAPDFDIR/lib/libLHAPDF.a";
          CONDITIONAL_LHAPDFINCS="-I$CONDITIONAL_LHAPDFDIR/include";
        fi;
        if test -f "$CONDITIONAL_LHAPDFDIR/include/LHAPDF/LHAPDF.h"; then
          lhapdfnativewrapper=true;
        fi;
        AC_MSG_RESULT([${CONDITIONAL_LHAPDFDIR}]); lhapdf=true;
      else
        AC_MSG_ERROR(Unable to determine path to your LHAPDF installation. \
                     Please set \$LHAPDFDIR manually.);
      fi;
    ],
    [ lhapdf=false ]
  )
  if test "$lhapdf" = "true" ; then
    AC_DEFINE_UNQUOTED([LHAPDF_PATH], "$CONDITIONAL_LHAPDFDIR", [LHAPDF directory])
    AC_DEFINE([USING__LHAPDF], "1", [using LHAPDF])
  fi
  AC_SUBST(CONDITIONAL_LHAPDFDIR)
  AC_SUBST(CONDITIONAL_LHAPDFLIBS)
  AC_SUBST(CONDITIONAL_LHAPDFINCS)
  AM_CONDITIONAL(LHAPDF_SUPPORT, test "$lhapdf" = "true")

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
                CONDITIONAL_CERNLIBLIBS="${enableval}/lib/libpacklib_noshift.a ${enableval}/lib/libmathlib.a ${enableval}/lib/libkernlib_noshift.a"
                cernlib=true;
                AC_MSG_RESULT(${enableval});
              elif test -f "${enableval}/lib/libkernlib.a"; then
                CONDITIONAL_CERNLIBLIBS="${enableval}/lib/libpacklib.a ${enableval}/lib/libmathlib.a ${enableval}/lib/libkernlib.a"
                cernlib=true;
                AC_MSG_RESULT(${enableval});
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
      esac ],
    [ zlib=false ]
  )
  if test "$zlib" = "true" ; then
    AC_DEFINE([USING__GZIP], "1", [using gzip])
  fi
  AM_CONDITIONAL(GZIP_SUPPORT, test "$zlib" = "true")
  AC_SUBST(CONDITIONAL_GZIPLIBS)

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

  AC_ARG_ENABLE(amisicinclude,
    AC_HELP_STRING([--disable-amisicinclude], [Disable inclusion of AMISIC headers]),
    [ AC_MSG_CHECKING(whether to include AMISIC headers);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); amisicinclude=true;;
        no)  AC_MSG_RESULT(no); amisicinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include AMISIC stuff); AC_MSG_RESULT(yes); amisicinclude=true; ]
  )
  if test "$amisicinclude" = "true" ; then
    AC_DEFINE([USING__Amisic], "1", [using AMISIC])
    CONDITIONAL_AMISICLIBS="\${AMISICLIBS}"
  fi
  AC_SUBST(CONDITIONAL_AMISICLIBS)
  AM_CONDITIONAL(AMISIC_SUPPORT, test "$amisicinclude" = "true" )
  
  AC_ARG_ENABLE(ahadicinclude,
    AC_HELP_STRING([--disable-ahadicinclude], [Disable inclusion of AHADIC headers]),
    [ AC_MSG_CHECKING(whether to include AHADIC headers);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); ahadicinclude=true;;
        no)  AC_MSG_RESULT(no); ahadicinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include AHADIC stuff); AC_MSG_RESULT(yes); ahadicinclude=true; ]
  )
  if test "$ahadicinclude" = "true" ; then
    AC_DEFINE([USING__Ahadic], "1", [using AHADIC])
    CONDITIONAL_AHADICLIBS="\${AHADICLIBS}"
  fi
  AC_SUBST(CONDITIONAL_AHADICLIBS)
  AM_CONDITIONAL(AHADIC_SUPPORT, test "$ahadicinclude" = "true" )
  
  AC_ARG_ENABLE(hadronsinclude,
    AC_HELP_STRING([--disable-hadronsinclude], [Disable inclusion of HADRONS headers]),
    [ AC_MSG_CHECKING(whether to include HADRONS headers);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); hadronsinclude=true;;
        no)  AC_MSG_RESULT(no); hadronsinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include HADRONS stuff); AC_MSG_RESULT(yes); hadronsinclude=true; ]
  )
  if test "$hadronsinclude" = "true" ; then
    AC_DEFINE([USING__Hadrons], "1", [using HADRONS])
    CONDITIONAL_HADRONSLIBS="\${HADRONSLIBS}"
  fi
  AC_SUBST(CONDITIONAL_HADRONSLIBS)
  AM_CONDITIONAL(HADRONS_SUPPORT, test "$hadronsinclude" = "true" )
  
  AC_ARG_ENABLE(photonsinclude,
    AC_HELP_STRING([--disable-photonsinclude], [Disable inclusion of PHOTONS headers]),
    [ AC_MSG_CHECKING(whether to include PHOTONS headers);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); photonsinclude=true;;
        no)  AC_MSG_RESULT(no); photonsinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include PHOTONS stuff); AC_MSG_RESULT(yes); photonsinclude=true; ]
  )
  if test "$photonsinclude" = "true" ; then
    AC_DEFINE([USING__Photons], "1", [using PHOTONS])
    CONDITIONAL_PHOTONSLIBS="\${PHOTONSLIBS}"
  fi
  AC_SUBST(CONDITIONAL_PHOTONSLIBS)
  AM_CONDITIONAL(PHOTONS_SUPPORT, test "$photonsinclude" = "true" )
])
