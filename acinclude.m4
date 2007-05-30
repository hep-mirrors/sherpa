dnl set flags according to build environment

AC_DEFUN([SHERPA_SETUP_BUILDSYSTEM],
[
  case "$build_os:$build_cpu:$build_vendor" in
    *darwin*:*power*:*)
      echo "checking for architecture... Darwin MacOS"
      ldflags="-dynamic -flat_namespace"
      AC_DEFINE([ARCH_DARWIN], "1", [Architecture identified as Darwin MacOS]) ;;
    *linux*:*:*)
      echo "checking for architecture...  Linux"
      ldflags="-rdynamic"
      AC_DEFINE([ARCH_LINUX], "1", [Architecture identified as Linux]) ;;
    *)
      echo "checking for architecture...  unknown"
      echo "hosts system type $build not yet supported, assuming unix behaviour."
      echo "possible failure due to unknown compiler/linker characteristics."
      echo "please inform us about build results at info@sherpa-mc.de"
      echo "(will continue in 10 seconds)"
      sleep 10
      ldflags="-rdynamic"
      AC_DEFINE([ARCH_UNIX], "1", [Architecture identified as Unix]) ;;
  esac
  AC_SUBST(ldflags)
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
  AMEGICDIR="\${top_srcdir}/AMEGIC++-2.0"
  AMEGICBUILDDIR="\${top_builddir}/AMEGIC++-2.0"
  AMEGICINCS="-I\${AMEGICDIR}/Main -I\${AMEGICDIR}/Amplitude -I\${AMEGICDIR}/Phasespace \
              -I\${AMEGICDIR}/String -I\${AMEGICDIR}/Model -I\${AMEGICDIR}/Amplitude/Zfunctions"
  AMEGICLIBS="-L\${AMEGICBUILDDIR}/Main -L\${AMEGICBUILDDIR}/Amplitude -L\${AMEGICBUILDDIR}/Phasespace \
              -L\${AMEGICBUILDDIR}/String -L\${AMEGICBUILDDIR}/Model -L\${AMEGICBUILDDIR}/Amplitude/Zfunctions \        
              -lAmegic -lAmplitude -lAmegicPSGen -lZfunctions -lModel -lString"
  AC_SUBST(AMEGICDIR)
  AC_SUBST(AMEGICBUILDDIR)
  AC_SUBST(AMEGICINCS)
  AC_SUBST(AMEGICLIBS)

  AMISICDIR="\${top_srcdir}/AMISIC++-1.0"
  AMISICBUILDDIR="\${top_builddir}/AMISIC++-1.0"
  AMISICINCS="-I\${AMISICDIR}/Main -I\${AMISICDIR}/Tools -I\${AMISICDIR}/Model"
  AMISICLIBS="-L\${AMISICBUILDDIR}/Main -L\${AMISICBUILDDIR}/Tools -L\${AMISICBUILDDIR}/Model \
              -lAmisic -lAmisicModel -lAmisicTools"
  AC_SUBST(AMISICDIR)
  AC_SUBST(AMISICBUILDDIR)
  AC_SUBST(AMISICINCS)
  AC_SUBST(AMISICLIBS)

  AHADICDIR="\${top_srcdir}/AHADIC++-1.0"
  AHADICBUILDDIR="\${top_builddir}/AHADIC++-1.0"
  AHADICINCS="-I\${AHADICDIR}/Main -I\${AHADICDIR}/Tools -I\${AHADICDIR}/Formation \
	      -I\${AHADICDIR}/Decays"
  AHADICLIBS="-L\${AHADICBUILDDIR}/Main -L\${AHADICBUILDDIR}/Tools -L\${AHADICBUILDDIR}/Formation -L\${AHADICBUILDDIR}/Decays \
              -lAhadicMain -lAhadicTools -lAhadicFormation -lAhadicDecays"
  AC_SUBST(AHADICDIR)
  AC_SUBST(AHADICBUILDDIR)
  AC_SUBST(AHADICINCS)
  AC_SUBST(AHADICLIBS)
  
  ANALYSISDIR="\${top_srcdir}/ANALYSIS-1.0"
  ANALYSISBUILDDIR="\${top_builddir}/ANALYSIS-1.0"
  ANALYSISINCS="-I\${ANALYSISDIR}/Tools -I\${ANALYSISDIR}/Main -I\${ANALYSISDIR}/Triggers -I\${ANALYSISDIR}/Observables"
  ANALYSISLIBS="-L\${ANALYSISBUILDDIR}/Tools -L\${ANALYSISBUILDDIR}/Main -L\${ANALYSISBUILDDIR}/Triggers -L\${ANALYSISBUILDDIR}/Observables -lAnalysis -lAnalysisTools -lAnalysisTriggers -lObservables"
  AC_SUBST(ANALYSISDIR)
  AC_SUBST(ANALYSISBUILDDIR)
  AC_SUBST(ANALYSISINCS)
  AC_SUBST(ANALYSISLIBS)
  
  APACICDIR="\${top_srcdir}/APACIC++-2.0"
  APACICBUILDDIR="\${top_builddir}/APACIC++-2.0"
  APACICINCS="-I\${APACICDIR}/Main -I\${APACICDIR}/Showers"
  APACICLIBS="-L\${APACICBUILDDIR}/Main -L\${APACICBUILDDIR}/Showers -lApacicShowers -lApacicMain"
  AC_SUBST(APACICDIR)
  AC_SUBST(APACICBUILDDIR)
  AC_SUBST(APACICINCS)
  AC_SUBST(APACICLIBS)
  
  ATOOLSDIR="\${top_srcdir}/ATOOLS-2.0"
  ATOOLSBUILDDIR="\${top_builddir}/ATOOLS-2.0"
  ATOOLSINCS="-I\${ATOOLSDIR}/Phys -I\${ATOOLSDIR}/Math -I\${ATOOLSDIR}/Org"
  ATOOLSLIBS="-L\${ATOOLSBUILDDIR}/Phys -L\${ATOOLSBUILDDIR}/Math -L\${ATOOLSBUILDDIR}/Org -lToolsPhys -lToolsMath -lToolsOrg"
  AC_SUBST(ATOOLSDIR)
  AC_SUBST(ATOOLSBUILDDIR)
  AC_SUBST(ATOOLSINCS)
  AC_SUBST(ATOOLSLIBS)
  
  BEAMDIR="\${top_srcdir}/BEAM-1.0"
  BEAMBUILDDIR="\${top_builddir}/BEAM-1.0"
  BEAMINCS="-I\${BEAMDIR}/Main"
  BEAMLIBS="-L\${BEAMBUILDDIR}/Main -lBeam"
  AC_SUBST(BEAMDIR)
  AC_SUBST(BEAMBUILDDIR)
  AC_SUBST(BEAMINCS)
  AC_SUBST(BEAMLIBS)
  
  EXTRAXSDIR="\${top_srcdir}/EXTRA_XS-1.0"
  EXTRAXSBUILDDIR="\${top_builddir}/EXTRA_XS-1.0"
  EXTRAXSINCS="-I\${EXTRAXSDIR}/COBG -I\${EXTRAXSDIR}/CDBG -I\${EXTRAXSDIR}/BFKL -I\${EXTRAXSDIR}/Two2Two -I\${EXTRAXSDIR}/Main"
  EXTRAXSLIBS="-L\${EXTRAXSBUILDDIR}/COBG -L\${EXTRAXSBUILDDIR}/CDBG -L\${EXTRAXSBUILDDIR}/BFKL -L\${EXTRAXSBUILDDIR}/Two2Two -L\${EXTRAXSBUILDDIR}/Main -lExtraXS -lExtraXS2_2 -lExtraXSBFKL -lExtraXSCOBG -lExtraXSCDBG"
  AC_SUBST(EXTRAXSDIR)
  AC_SUBST(EXTRAXSBUILDDIR)
  AC_SUBST(EXTRAXSINCS)
  AC_SUBST(EXTRAXSLIBS)
  
  HADRONSDIR="\${top_srcdir}/HADRONS++-0.0"
  HADRONSBUILDDIR="\${top_builddir}/HADRONS++-0.0"
  HADRONSINCS="-I\${HADRONSDIR}/Main -I\${HADRONSDIR}/ME_Library -I\${HADRONSDIR}/PS_Library"
  HADRONSLIBS="-L\${HADRONSBUILDDIR}/Main -L\${HADRONSBUILDDIR}/ME_Library -L\${HADRONSBUILDDIR}/PS_Library \
               -lHadronsMain -lHadronsMEs -lHadronsPSs"
  AC_SUBST(HADRONSDIR)
  AC_SUBST(HADRONSBUILDDIR)
  AC_SUBST(HADRONSINCS)
  AC_SUBST(HADRONSLIBS)
  
  MODELDIR="\${top_srcdir}/MODEL-1.0"
  MODELBUILDDIR="\${top_builddir}/MODEL-1.0"
  MODELINCS="-I\${MODELDIR}/Main"
  MODELLIBS="-L\${MODELBUILDDIR}/Main -lModelMain"
  AC_SUBST(MODELDIR)
  AC_SUBST(MODELBUILDDIR)
  AC_SUBST(MODELINCS)
  AC_SUBST(MODELLIBS)
  
  HDECAYINCS="-I\${MODELDIR}/Hdecay"
  HDECAYLIBS="-L\${MODELBUILDDIR}/Hdecay -lHdecay"
  AC_SUBST(HDECAYINCS)
  AC_SUBST(HDECAYLIBS)
  
  ISAJETINCS="-I\${MODELDIR}/Isajet"
  ISAJETLIBS="-L\${MODELBUILDDIR}/Isajet -lIsajet"
  AC_SUBST(ISAJETINCS)
  AC_SUBST(ISAJETLIBS)
  
  PDFDIR="\${top_srcdir}/PDF-1.0"
  PDFBUILDDIR="\${top_builddir}/PDF-1.0"
  PDFINCS="-I\${PDFDIR}/Main -I\${PDFDIR}/Remnant -I\${PDFDIR}/LHAPDF -I\${PDFDIR}/MRST \
           -I\${PDFDIR}/GRV -I\${PDFDIR}/Sudakov -I\${PDFDIR}/KMR -I\${PDFDIR}/Remnant"
  PDFLIBS="-L\${PDFBUILDDIR}/Main -L\${PDFBUILDDIR}/Remnant -L\${PDFBUILDDIR}/LHAPDF -L\${PDFBUILDDIR}/MRST \
           -L\${PDFBUILDDIR}/GRV -L\${PDFBUILDDIR}/Sudakov -L\${PDFBUILDDIR}/KMR -L\${PDFBUILDDIR}/Remnant \
           -lPDF -lDUPDF -lSudakov -lMRST -lGRV -lRemnant"
  AC_SUBST(PDFDIR)
  AC_SUBST(PDFBUILDDIR)
  AC_SUBST(PDFINCS)
  AC_SUBST(PDFLIBS)
  
  CTEQINCS="-I\${PDFDIR}/CTEQ"
  CTEQLIBS="-L\${PDFBUILDDIR}/CTEQ -lCTEQ"
  AC_SUBST(CTEQINCS)
  AC_SUBST(CTEQLIBS)
  
  PHASICDIR="\${top_srcdir}/PHASIC++-1.0"
  PHASICBUILDDIR="\${top_builddir}/PHASIC++-1.0"
  PHASICINCS="-I\${PHASICDIR}/Main -I\${PHASICDIR}/Foam"
  PHASICLIBS="-L\${PHASICBUILDDIR}/Main -lPhasespace"
  AC_SUBST(PHASICDIR)
  AC_SUBST(PHASICBUILDDIR)
  AC_SUBST(PHASICINCS)
  AC_SUBST(PHASICLIBS)
  
  SHERPADIR="\${top_srcdir}/SHERPA-1.0"
  SHERPABUILDDIR="\${top_builddir}/SHERPA-1.0"
  SHERPAINCS="-I\${SHERPADIR}/Single_Events -I\${SHERPADIR}/PerturbativePhysics \
              -I\${SHERPADIR}/LundTools -I\${SHERPADIR}/Tools -I\${SHERPADIR}/Main \
              -I\${SHERPADIR}/Initialization -I\${SHERPADIR}/SoftPhysics -I\${SHERPADIR}/HerwigTools"
  SHERPALIBS="-L\${SHERPABUILDDIR}/Single_Events -L\${SHERPABUILDDIR}/PerturbativePhysics \
              -L\${SHERPABUILDDIR}/LundTools -L\${SHERPABUILDDIR}/Tools -L\${SHERPABUILDDIR}/Main \
              -L\${SHERPABUILDDIR}/Initialization -L\${SHERPABUILDDIR}/SoftPhysics -L\${SHERPABUILDDIR}/HerwigTools \
              -lSherpaMain -lSherpaInitialization -lSherpaSingleEvents \
              -lSherpaPerturbativePhysics -lSherpaSoftPhysics -lLundTools -lSherpaTools"
  SHERPAFLAGS="-pedantic -Wall"
  AC_SUBST(SHERPADIR)
  AC_SUBST(SHERPABUILDDIR)
  AC_SUBST(SHERPAINCS)
  AC_SUBST(SHERPALIBS)
  AC_SUBST(SHERPAFLAGS)

  if test "x$prefix" = "xNONE"; then
    prefix=$ac_default_prefix
  fi
  if test "x$exec_prefix" = "xNONE"; then
    exec_prefix=$prefix
  fi

  AS_AC_EXPAND(LIBDIR, ${libdir})
  AS_AC_EXPAND(INCLUDEDIR, ${includedir})
  AS_AC_EXPAND(BINDIR, ${bindir})
  AS_AC_EXPAND(DATADIR, ${datadir})

  AC_DEFINE_UNQUOTED([SHERPA_VERSION], ["`echo AC_PACKAGE_VERSION | cut -d. -f1`"], [Sherpa version])
  AC_DEFINE_UNQUOTED([SHERPA_SUBVERSION], ["`echo AC_PACKAGE_VERSION | cut -d. -f2,3`"], [Sherpa subversion])
  AC_DEFINE_UNQUOTED([SHERPA_BUILD_PATH], "$PWD", [Sherpa build path])
  AC_DEFINE_UNQUOTED([SHERPA_INCLUDE_PATH], "$INCLUDEDIR/SHERPA-MC", [Sherpa include directory])
  AC_DEFINE_UNQUOTED([SHERPA_BINARY_PATH], "$BINDIR/SHERPA-MC", [Sherpa binary directory])
  AC_DEFINE_UNQUOTED([SHERPA_LIBRARY_PATH], "$LIBDIR/SHERPA-MC", [Sherpa library directory])
  AC_DEFINE_UNQUOTED([SHERPA_PDFS_PATH], "$DATADIR/SHERPA-MC", [Sherpa data directory])
  AC_DEFINE_UNQUOTED([SHERPA_SHARE_PATH], "$DATADIR/SHERPA-MC", [Sherpa data directory])
  AC_DEFINE([USING__COLOUR], "1", [Using colour])
])



dnl Conditional compiling and linking

AC_DEFUN([SHERPA_SETUP_CONFIGURE_OPTIONS],
[
  AC_ARG_ENABLE(
    svninclude,
    AC_HELP_STRING([--disable-svninclude], [Don't distribute SVN synchronization directories.]),
    [ AC_MSG_CHECKING(whether to enable SVN synchronization)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no);
              SVNINCLUDE="";;
        yes) AC_MSG_RESULT(yes);
              SVNINCLUDE=".svn";;
      esac ],
    [ AC_MSG_CHECKING(whether to enable SVN synchronization); AC_MSG_RESULT(yes); SVNINCLUDE=".svn" ] 
  )
  AC_SUBST(SVNINCLUDE)

  AC_ARG_ENABLE(
    mcatnloinclude,
    AC_HELP_STRING([--enable-mcatnloinclude], [Enable MC@NLO support]),
    [ AC_MSG_CHECKING(for MC@NLO support)
      case "${enableval}" in
        no)  AC_MSG_RESULT(no); mcatnloinclude=false ;;
        yes) AC_MSG_RESULT(yes); mcatnloinclude=true ;;
      esac ],
    [ AC_MSG_CHECKING(for MC@NLO support); AC_MSG_RESULT(no); mcatnloinclude=false ] 
  )
  if test "$mcatnloinclude" = "true" ; then
    AC_DEFINE([USING__MCatNLO], "1", [using MC@NLO])
    CONDITIONAL_MCATNLOLIBS="\${MCATNLOLIBS}"
  fi
  AC_SUBST(CONDITIONAL_MCATNLOLIBS)
  AM_CONDITIONAL(MCATNLO_SUPPORT, test "$mcatnloinclude" = "true" )
  
  AC_ARG_ENABLE(
    clhep,
    AC_HELP_STRING([--enable-clhep], [Enable CLHEP support]),
    [ AC_MSG_CHECKING(for CLHEP installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(CLHEP not enabled); clhep=false ;;
        yes) if test -d "$CLHEPDIR"; then
                CONDITIONAL_CLHEPDIR="$CLHEPDIR"
                CONDITIONAL_CLHEPINCS="-I$CLHEPDIR/include"
                possible_libs="libCLHEP-g++.*.a libCLHEP-g++.*.so libCLHEP.so libCLHEP-1*.so libCLHEP-2*.so"
                for J in $possible_libs; do
                  result=`find $CLHEPDIR/lib -name "$J" | head -n 1`;
                  if test "$result" != ""; then
                    result=`basename $result | sed -e 's/lib//' | sed -e 's/\.so//' | sed -e 's/\.a//'`
                    break;
                  fi
                done;
                if test "$result" != ""; then
                  CONDITIONAL_CLHEPLIBS="-L$CLHEPDIR/lib -R$CLHEPDIR/lib -l$result";
                else
                  AC_MSG_ERROR(Did not find any library of the following type in $CLHEPDIR/lib: $possible_libs and clhep-config was not available in \$PATH.);
                fi
              elif test -x "`which clhep-config`"; then
                CONDITIONAL_CLHEPDIR=`clhep-config --prefix`;
                CONDITIONAL_CLHEPLIBS=`clhep-config --libs`
                CONDITIONAL_CLHEPINCS=`clhep-config --include`
                if ! test -d "$CONDITIONAL_CLHEPDIR"; then
                  AC_MSG_ERROR(clhep-config --prefix returned a path that is not available. Please check your CLHEP installation and set \$CLHEPDIR manually.);
                fi
              else
                AC_MSG_ERROR(\$CLHEPDIR is not a valid path and clhep-config was not found.);
              fi;
              AC_MSG_RESULT([${CONDITIONAL_CLHEPDIR}]); clhep=true;;
      esac ],
    [ clhep=false ]
  )
  if test "$clhep" = "true" ; then
    AC_DEFINE([USING__CLHEP], "1", [Using CLHEP])
  fi
  AC_SUBST(CONDITIONAL_CLHEPDIR)
  AC_SUBST(CONDITIONAL_CLHEPINCS)
  AC_SUBST(CONDITIONAL_CLHEPLIBS)
  AM_CONDITIONAL(CLHEP_SUPPORT, test "$clhep" = "true")

  AC_ARG_ENABLE(
    hepmc2,
    AC_HELP_STRING([--enable-hepmc2], [Enable HepMC (version 2.x) support]),
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
      esac ],
    [ hepmc2=false ]
  )
  if test "$hepmc2" = "true" ; then
    AC_DEFINE([USING__HEPMC2], "1", [Using HEPMC2])
  fi
  AC_SUBST(CONDITIONAL_HEPMC2DIR)
  AC_SUBST(CONDITIONAL_HEPMC2INCS)
  AC_SUBST(CONDITIONAL_HEPMC2LIBS)
  AM_CONDITIONAL(HEPMC2_SUPPORT, test "$hepmc2" = "true")

  AC_ARG_ENABLE(
    root,
    AC_HELP_STRING([--enable-root], [Enable ROOT support]),
    [ AC_MSG_CHECKING(for ROOT installation directory)
      case "${enableval}" in
        no)  AC_MSG_RESULT(ROOT not enabled); root=false;;
        yes) if test -d "$ROOTSYS"; then
               CONDITIONAL_ROOTDIR=$ROOTSYS
               CONDITIONAL_ROOTINCS=-I`$ROOTSYS/bin/root-config --incdir`;
               CONDITIONAL_ROOTLIBS=`$ROOTSYS/bin/root-config --glibs`
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
    AC_HELP_STRING([--enable-lhapdf], [Enable LHAPDF support]),
    [ AC_MSG_CHECKING(for LHAPDF installation directory);
      case "${enableval}" in
        no)  AC_MSG_RESULT(LHAPDF not enabled); lhapdf=false ;;
        yes) if test -d "$LHAPDFDIR"; then
              CONDITIONAL_LHAPDFDIR=$LHAPDFDIR;
              CONDITIONAL_LHAPDFLIBS="-lLHAPDF $LHAPDFDIR/lib/libLHAPDF.a"
            elif test -x "`which lhapdf-config`"; then
              CONDITIONAL_LHAPDFDIR=`lhapdf-config --prefix`;
              CONDITIONAL_LHAPDFLIBS="-lLHAPDF `lhapdf-config --prefix`/lib/libLHAPDF.a"
              if ! test -d "$CONDITIONAL_LHAPDFDIR"; then
                AC_MSG_ERROR(lhapdf-config --prefix returned a path that is not available. Please check your LHAPDF installation and set \$LHAPDFDIR manually.);
              fi
  else
              AC_MSG_ERROR(\$LHAPDFDIR is not a valid path and lhapdf-config was not found.);
            fi;
            AC_MSG_RESULT([${CONDITIONAL_LHAPDFDIR}]); lhapdf=true;;
      esac ],
    [ lhapdf=false ]
  )
  if test "$lhapdf" = "true" ; then
    AC_DEFINE([USING__LHAPDF], "1", [using LHAPDF])
  fi
  AC_SUBST(CONDITIONAL_LHAPDFDIR)
  AC_SUBST(CONDITIONAL_LHAPDFLIBS)
  AM_CONDITIONAL(LHAPDF_SUPPORT, test "$lhapdf" = "true")

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

  AC_ARG_ENABLE(modelinclude,
    AC_HELP_STRING([--disable-modelinclude], [Disable inclusion of MODEL headers]),
    [ AC_MSG_CHECKING(whether to include MODEL headers);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); modelinclude=true;;
        no)  AC_MSG_RESULT(no); modelinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include MODEL headers); AC_MSG_RESULT(yes); modelinclude=true; ]
  )
  if test "$modelinclude" = "true" ; then
    AC_DEFINE([USING__Model], "1", [using Model])
  else
    AC_DEFINE([USING__ATOOLS_only], "1", [not using Model, using ATOOLS_only])
  fi
  AM_CONDITIONAL(MODEL_SUPPORT, test "$modelinclude" = "true" )

  AC_ARG_ENABLE(isajetinclude,
    AC_HELP_STRING([--disable-isajetinclude], [Disable inclusion of Isajet stuff]),
    [ AC_MSG_CHECKING(whether to include Isajet stuff);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); isajetinclude=true;;
        no)  AC_MSG_RESULT(no); isajetinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include Isajet stuff); AC_MSG_RESULT(yes); isajetinclude=true; ]
  )
  dnl comment out if isajet to be compiled
  echo "Hardwired: Omitting all Isajet stuff"; isajetinclude=false;
  if test "$isajetinclude" = "true" ; then
    AC_DEFINE([USING__Isajet], "1", [using Isajet])
    CONDITIONAL_ISAJETINCS="\${ISAJETINCS}"
    CONDITIONAL_ISAJETLIBS="\${ISAJETLIBS}"
  fi
  AC_SUBST(CONDITIONAL_ISAJETINCS)
  AC_SUBST(CONDITIONAL_ISAJETLIBS)
  AM_CONDITIONAL(ISAJET_SUPPORT, test "$isajetinclude" = "true" )
  
  AC_ARG_ENABLE(hdecayinclude,
    AC_HELP_STRING([--disable-hdecayinclude], [Disable inclusion of Hdecay stuff]),
    [ AC_MSG_CHECKING(whether to include Hdecay stuff);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); hdecayinclude=true;;
        no)  AC_MSG_RESULT(no); hdecayinclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include Hdecay stuff); AC_MSG_RESULT(yes); hdecayinclude=true; ]
  )
  dnl comment out if hdecay to be compiled
  echo "Hardwired: Omitting all Hdecay stuff"; hdecayinclude=false;
  if test "$hdecayinclude" = "true" ; then
    AC_DEFINE([USING__Hdecay], "1", [using Hdecay])
    CONDITIONAL_HDECAYLIBS="\${HDECAYLIBS}"
    CONDITIONAL_HDECAYINCS="\${HDECAYINCS}"
  fi
  AC_SUBST(CONDITIONAL_HDECAYLIBS)
  AC_SUBST(CONDITIONAL_HDECAYINCS)
  AM_CONDITIONAL(HDECAY_SUPPORT, test "$hdecayinclude" = "true" )
  
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
    CONDITIONAL_AMISICINCS="\${AMISICINCS}"
  fi
  AC_SUBST(CONDITIONAL_AMISICLIBS)
  AC_SUBST(CONDITIONAL_AMISICINCS)
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
    CONDITIONAL_AHADICINCS="\${AHADICINCS}"
  fi
  AC_SUBST(CONDITIONAL_AHADICLIBS)
  AC_SUBST(CONDITIONAL_AHADICINCS)
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
    CONDITIONAL_HADRONSINCS="\${HADRONSINCS}"
  fi
  AC_SUBST(CONDITIONAL_HADRONSLIBS)
  AC_SUBST(CONDITIONAL_HADRONSINCS)
  AM_CONDITIONAL(HADRONS_SUPPORT, test "$hadronsinclude" = "true" )
  
  SHERPASUBDIR="SHERPA-1.0"
  AC_ARG_ENABLE(sherpainclude,
    AC_HELP_STRING([--disable-sherpainclude], [Disable inclusion of SHERPA headers]),
    [ AC_MSG_CHECKING(whether to include SHERPA headers);
      case "${enableval}" in
        yes) AC_MSG_RESULT(yes); sherpainclude=true;;
        no)  AC_MSG_RESULT(no); sherpainclude=false;;
      esac ],
    [ AC_MSG_CHECKING(whether to include SHERPA stuff); AC_MSG_RESULT(yes); sherpainclude=true; ]
  )
  if test "$sherpainclude" = "true" ; then
    AC_DEFINE([USING__Sherpa], "1", [using SHERPA])
    CONDITIONAL_SHERPAINCS="\${SHERPAINCS}"
    CONDITIONAL_SHERPALIBS="\${SHERPALIBS}"
  fi
  AC_SUBST(CONDITIONAL_SHERPALIBS)
  AC_SUBST(CONDITIONAL_SHERPAINCS)
  AM_CONDITIONAL(SHERPA_SUPPORT, test "$sherpainclude" = "true" )
])
