dnl workaround for old automake on darwin

AC_DEFUN([AM_CONFIG_HEADERS], [AC_CONFIG_HEADERS($@)])

dnl set flags according to build environment

AC_DEFUN([SHERPA_SETUP_BUILDSYSTEM],
[
  case "$build_os:$build_cpu:$build_vendor" in
    *darwin*:*:*)
      echo "checking for architecture... Darwin MacOS"
      ldflags="-dynamic -flat_namespace"
      AC_DEFINE([ARCH_DARWIN], "1", [Architecture identified as Darwin MacOS])
      AC_DEFINE([LIB_SUFFIX], ".dylib", [library suffix set to .dylib]) 
      AC_DEFINE([LD_PATH_NAME], "DYLD_LIBRARY_PATH", [ld path name set to DYLD_LIBRARY_PATH]) ;;
    *linux*:*:*)
      echo "checking for architecture...  Linux"
      ldflags="-rdynamic"
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
      ldflags="-rdynamic"
      AC_DEFINE([ARCH_UNIX], "1", [Architecture identified as Unix])
      AC_DEFINE([LIB_SUFFIX], ".so", [library suffix set to .so]) 
      AC_DEFINE([LD_PATH_NAME], "LD_LIBRARY_PATH", [ld path name set to LD_LIBRARY_PATH]) ;;
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
  case "$build_os:$build_cpu:$build_vendor" in
    *darwin*:*:*)
      CFL=$(echo $FLIBS | awk '{ for (i=1;i<NF;++i) \
        if (match($i,"-lSystem")==0 && match($i,"-lgcc_s")==0) printf " "$i; }')
      FLIBS=$CFL
      echo "trimming fortran libs for Darwin... "$FLIBS
      AC_SUBST(FLIBS)
      f77_main_darwin=`(test $F77 != g77 && echo MAIN__) || echo main`
      echo "setting fortran main name for Darwin... "$f77_main_darwin
      AC_DEFINE_UNQUOTED([F77_MAIN], [`echo $f77_main_darwin`],[alternate entry point]) ;;
    *)
      ;;
  esac

  AMEGICDIR="\${top_srcdir}/AMEGIC++"
  AMEGICBUILDDIR="\${top_builddir}/AMEGIC++"
  AMEGICINCS="-I\${AMEGICDIR}/Main -I\${AMEGICDIR}/Amplitude -I\${AMEGICDIR}/Phasespace \
              -I\${AMEGICDIR}/String -I\${AMEGICDIR}/Amplitude/Zfunctions"
  AMEGICLIBS="-L\${AMEGICBUILDDIR}/Main -L\${AMEGICBUILDDIR}/Amplitude -L\${AMEGICBUILDDIR}/Phasespace \
              -L\${AMEGICBUILDDIR}/String -L\${AMEGICBUILDDIR}/Amplitude/Zfunctions \        
              -lAmegic -lAmplitude -lAmegicPSGen -lZfunctions -lString"
  AC_SUBST(AMEGICDIR)
  AC_SUBST(AMEGICBUILDDIR)
  AC_SUBST(AMEGICINCS)
  AC_SUBST(AMEGICLIBS)

  AMISICDIR="\${top_srcdir}/AMISIC++"
  AMISICBUILDDIR="\${top_builddir}/AMISIC++"
  AMISICINCS="-I\${AMISICDIR}/Main -I\${AMISICDIR}/Tools -I\${AMISICDIR}/Model"
  AMISICLIBS="-L\${AMISICBUILDDIR}/Main -L\${AMISICBUILDDIR}/Tools -L\${AMISICBUILDDIR}/Model \
              -lAmisic -lAmisicModel -lAmisicTools"
  AC_SUBST(AMISICDIR)
  AC_SUBST(AMISICBUILDDIR)
  AC_SUBST(AMISICINCS)
  AC_SUBST(AMISICLIBS)

  AHADICDIR="\${top_srcdir}/AHADIC++"
  AHADICBUILDDIR="\${top_builddir}/AHADIC++"
  AHADICINCS="-I\${AHADICDIR}/Main -I\${AHADICDIR}/Tools -I\${AHADICDIR}/Formation \
	      -I\${AHADICDIR}/Decays"
  AHADICLIBS="-L\${AHADICBUILDDIR}/Main -L\${AHADICBUILDDIR}/Tools -L\${AHADICBUILDDIR}/Formation -L\${AHADICBUILDDIR}/Decays \
              -lAhadicMain -lAhadicTools -lAhadicFormation -lAhadicDecays"
  AC_SUBST(AHADICDIR)
  AC_SUBST(AHADICBUILDDIR)
  AC_SUBST(AHADICINCS)
  AC_SUBST(AHADICLIBS)
  
  ANALYSISDIR="\${top_srcdir}/ANALYSIS"
  ANALYSISBUILDDIR="\${top_builddir}/ANALYSIS"
  ANALYSISINCS="-I\${ANALYSISDIR}/Tools -I\${ANALYSISDIR}/Detector -I\${ANALYSISDIR}/Main -I\${ANALYSISDIR}/Triggers -I\${ANALYSISDIR}/Observables"
  ANALYSISLIBS="-L\${ANALYSISBUILDDIR}/Tools -L\${ANALYSISBUILDDIR}/Main -L\${ANALYSISBUILDDIR}/Triggers -L\${ANALYSISBUILDDIR}/Detector -L\${ANALYSISBUILDDIR}/Observables -lAnalysis -lAnalysisTools -lAnalysisTriggers -lAnalysisDetector -lObservables"
  AC_SUBST(ANALYSISDIR)
  AC_SUBST(ANALYSISBUILDDIR)
  AC_SUBST(ANALYSISINCS)
  AC_SUBST(ANALYSISLIBS)
  
  APACICDIR="\${top_srcdir}/APACIC++"
  APACICBUILDDIR="\${top_builddir}/APACIC++"
  APACICINCS="-I\${APACICDIR}/Main -I\${APACICDIR}/Showers"
  APACICLIBS="-L\${APACICBUILDDIR}/Main -L\${APACICBUILDDIR}/Showers -lApacicShowers -lApacicMain"
  AC_SUBST(APACICDIR)
  AC_SUBST(APACICBUILDDIR)
  AC_SUBST(APACICINCS)
  AC_SUBST(APACICLIBS)
  
  ATOOLSDIR="\${top_srcdir}/ATOOLS"
  ATOOLSBUILDDIR="\${top_builddir}/ATOOLS"
  ATOOLSINCS="-I\${ATOOLSDIR}/Phys -I\${ATOOLSDIR}/Math -I\${ATOOLSDIR}/Org"
  ATOOLSLIBS="-L\${ATOOLSBUILDDIR}/Phys -L\${ATOOLSBUILDDIR}/Math \
              -L\${ATOOLSBUILDDIR}/Org \
              -lToolsPhys -lToolsMath -lToolsOrg"
  AC_SUBST(ATOOLSDIR)
  AC_SUBST(ATOOLSBUILDDIR)
  AC_SUBST(ATOOLSINCS)
  AC_SUBST(ATOOLSLIBS)
  
  BEAMDIR="\${top_srcdir}/BEAM"
  BEAMBUILDDIR="\${top_builddir}/BEAM"
  BEAMINCS="-I\${BEAMDIR}/Main"
  BEAMLIBS="-L\${BEAMBUILDDIR}/Main -lBeam"
  AC_SUBST(BEAMDIR)
  AC_SUBST(BEAMBUILDDIR)
  AC_SUBST(BEAMINCS)
  AC_SUBST(BEAMLIBS)

  HELICITIESDIR="\${top_srcdir}/HELICITIES"
  HELICITIESBUILDDIR="\${top_builddir}/HELICITIES"
  HELICITIESINCS="-I\${HELICITIESDIR}/Main"
  HELICITIESLIBS="-L\${HELICITIESBUILDDIR}/Main -lHelicitiesMain"
  AC_SUBST(HELICITIESDIR)
  AC_SUBST(HELICITIESBUILDDIR)
  AC_SUBST(HELICITIESINCS)
  AC_SUBST(HELICITIESLIBS)
  
  EXTRAXSDIR="\${top_srcdir}/EXTRA_XS"
  EXTRAXSBUILDDIR="\${top_builddir}/EXTRA_XS"
  EXTRAXSINCS="-I\${EXTRAXSDIR}/Two2Two -I\${EXTRAXSDIR}/Main -I\${EXTRAXSDIR}/Model"
  EXTRAXSLIBS="-L\${EXTRAXSBUILDDIR}/Two2Two -L\${EXTRAXSBUILDDIR}/Main -L\${EXTRAXSBUILDDIR}/Model -lExtraXSModel -lExtraXS -lExtraXS2_2"
  AC_SUBST(EXTRAXSDIR)
  AC_SUBST(EXTRAXSBUILDDIR)
  AC_SUBST(EXTRAXSINCS)
  AC_SUBST(EXTRAXSLIBS)
  
  HADRONSDIR="\${top_srcdir}/HADRONS++"
  HADRONSBUILDDIR="\${top_builddir}/HADRONS++"
  HADRONSINCS="-I\${HADRONSDIR}/Main -I\${HADRONSDIR}/ME_Library \
               -I\${HADRONSDIR}/Current_Library -I\${HADRONSDIR}/PS_Library"
  HADRONSLIBS="-L\${HADRONSBUILDDIR}/Main -L\${HADRONSBUILDDIR}/ME_Library \
               -L\${HADRONSBUILDDIR}/Current_Library -L\${HADRONSBUILDDIR}/PS_Library \
               -lHadronsMain -lHadronsMEs -lHadronsCurrents -lHadronsPSs"
  AC_SUBST(HADRONSDIR)
  AC_SUBST(HADRONSBUILDDIR)
  AC_SUBST(HADRONSINCS)
  AC_SUBST(HADRONSLIBS)
  
  PHOTONSDIR="\${top_srcdir}/PHOTONS++"
  PHOTONSBUILDDIR="\${top_builddir}/PHOTONS++"
  PHOTONSINCS="-I\${PHOTONSDIR}/Main -I\${PHOTONSDIR}/Tools -I\${PHOTONSDIR}/PhaseSpace \
               -I\${PHOTONSDIR}/MEs"
  PHOTONSLIBS="-L\${PHOTONSBUILDDIR}/Main -L\${PHOTONSBUILDDIR}/Tools \
               -L\${PHOTONSBUILDDIR}/PhaseSpace -L\${PHOTONSBUILDDIR}/MEs \
               -lPhotonsMain -lPhotonsTools -lPhotonsPhaseSpace -lPhotonsMEs"
  AC_SUBST(PHOTONSDIR)
  AC_SUBST(PHOTONSBUILDDIR)
  AC_SUBST(PHOTONSINCS)
  AC_SUBST(PHOTONSLIBS)
  
  MODELDIR="\${top_srcdir}/MODEL"
  MODELBUILDDIR="\${top_builddir}/MODEL"
  MODELINCS="-I\${MODELDIR}/Main -I\${MODELDIR}/Interaction_Models -I\${MODELDIR}/Decays"
  MODELLIBS="-L\${MODELBUILDDIR}/Main -L\${MODELBUILDDIR}/Interaction_Models \
	     -L\${MODELBUILDDIR}/Decays \	
             -lModelMain -lModelInteractions -lModelDecays"
  AC_SUBST(MODELDIR)
  AC_SUBST(MODELBUILDDIR)
  AC_SUBST(MODELINCS)
  AC_SUBST(MODELLIBS)
  
  PDFDIR="\${top_srcdir}/PDF"
  PDFBUILDDIR="\${top_builddir}/PDF"
  PDFINCS="-I\${PDFDIR}/Main -I\${PDFDIR}/Remnant -I\${PDFDIR}/LHAPDF -I\${PDFDIR}/MRST \
           -I\${PDFDIR}/GRV -I\${PDFDIR}/Sudakov -I\${PDFDIR}/Remnant"
  PDFLIBS="-L\${PDFBUILDDIR}/Main -L\${PDFBUILDDIR}/Remnant -L\${PDFBUILDDIR}/LHAPDF \
           -L\${PDFBUILDDIR}/GRV -L\${PDFBUILDDIR}/Sudakov -L\${PDFBUILDDIR}/Remnant \
           -lPDF -lSudakov -lGRV -lRemnant"
  AC_SUBST(PDFDIR)
  AC_SUBST(PDFBUILDDIR)
  AC_SUBST(PDFINCS)
  AC_SUBST(PDFLIBS)
  
  CTEQINCS="-I\${PDFDIR}/CTEQ"
  CTEQLIBS="-L\${PDFBUILDDIR}/CTEQ -lCTEQ"
  AC_SUBST(CTEQINCS)
  AC_SUBST(CTEQLIBS)

  MRSTLIBS="-L\${PDFBUILDDIR}/MRST -lMRST"
  AC_SUBST(MRSTLIBS)
  
  PHASICDIR="\${top_srcdir}/PHASIC++"
  PHASICBUILDDIR="\${top_builddir}/PHASIC++"
  PHASICINCS="-I\${PHASICDIR}/Main"
  PHASICLIBS="-L\${PHASICBUILDDIR}/Main -lPhasespace"
  AC_SUBST(PHASICDIR)
  AC_SUBST(PHASICBUILDDIR)
  AC_SUBST(PHASICINCS)
  AC_SUBST(PHASICLIBS)
  
  SHERPADIR="\${top_srcdir}/SHERPA"
  SHERPABUILDDIR="\${top_builddir}/SHERPA"
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
  AC_DEFINE_UNQUOTED([SHERPA_INCLUDE_PATH], "$INCLUDEDIR/SHERPA-MC", [Sherpa include directory])
  AC_DEFINE_UNQUOTED([SHERPA_LIBRARY_PATH], "$LIBDIR/SHERPA-MC", [Sherpa library directory])
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
  fi
  AC_SUBST(CONDITIONAL_HEPMC2DIR)
  AC_SUBST(CONDITIONAL_HEPMC2INCS)
  AC_SUBST(CONDITIONAL_HEPMC2LIBS)
  AM_CONDITIONAL(HEPMC2_SUPPORT, test "$hepmc2" = "true")

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
        else
          CONDITIONAL_LHAPDFLIBS="-L\${PDFBUILDDIR}/LHAPDF -lLHAPDFSherpa $CONDITIONAL_LHAPDFDIR/lib/libLHAPDF.a";
        fi;
        if test -f "$CONDITIONAL_LHAPDFDIR/include/LHAPDF/LHAPDF.h"; then
          lhapdfnativewrapper=true;
          CONDITIONAL_LHAPDFINCS="-I$CONDITIONAL_LHAPDFDIR/include";
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
    AC_DEFINE([USING__LHAPDF], "1", [using LHAPDF])
    if test "$lhapdfnativewrapper" = "true"; then
      AC_DEFINE([LHAPDF__NATIVE__WRAPPER], "1", [using native C++ wrapper])
    fi
  fi
  AC_SUBST(CONDITIONAL_LHAPDFDIR)
  AC_SUBST(CONDITIONAL_LHAPDFLIBS)
  AC_SUBST(CONDITIONAL_LHAPDFINCS)
  AM_CONDITIONAL(LHAPDF_SUPPORT, test "$lhapdf" = "true")
  AM_CONDITIONAL(LHAPDF_NATIVE_WRAPPER, test "$lhapdfnativewrapper" = "true")

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
    CONDITIONAL_PHOTONSINCS="\${PHOTONSINCS}"
  fi
  AC_SUBST(CONDITIONAL_PHOTONSLIBS)
  AC_SUBST(CONDITIONAL_PHOTONSINCS)
  AM_CONDITIONAL(PHOTONS_SUPPORT, test "$photonsinclude" = "true" )
])
