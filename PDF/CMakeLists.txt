
########################################################################
#
#  Automatically or semiautomaticall generated, do not edit.
#
########################################################################
# The following input was used
# CT14/Makefile.am
# Electron/Makefile.am
# GRS/Makefile.am
# GRV/Makefile.am
# Main/Makefile.am
# Makefile.am
# NNPDF/Makefile.am
# SAL/Makefile.am
# SASG/Makefile.am

########################################################################
#PDF/CT14/Makefile.am
 
# pkglib_LTLIBRARIES = libCT14Sherpa.la
#GITTAG = PDF/CT14 include ../../git.make libCT14Sherpa_la
  set(PDF_CT14_SOURCES  CT14Pdf.f CT14_Fortran_Interface.C)
#nobase_dist_pkgdata_DATA = CT14Grid/ct14n.00.pds CT14Grid/ct14nn.00.pds
  list(TRANSFORM PDF_CT14_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/CT14/")
 
#PDF/Electron/Makefile.am
 
# pkglib_LTLIBRARIES = libPDFESherpa.la
#GITTAG = PDF/Electron include ../../git.make libPDFESherpa_la
  set(PDF_Electron_SOURCES  PDF_Electron.C)
#localinc_HEADERS = PDF_Electron.H
  list(TRANSFORM PDF_Electron_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/Electron/")
 
#PDF/GRS/Makefile.am
 
# pkglib_LTLIBRARIES = libGRSSherpa.la
#GITTAG = PDF/GRS include ../../git.make AM_FFLAGS = -std = legacy libGRSSherpa_la
  set(PDF_GRS_SOURCES  grsg99.f GRSph_Fortran_Interface.C)
#localinc_HEADERS = GRSph_Fortran_Interface.H
#nobase_dist_pkgdata_DATA = GRSGrid/grsg99d.grid GRSGrid/grsg99lo.grid GRSGrid/grsg99m.grid GRSGrid/grsg99nlo.grid
  list(TRANSFORM PDF_GRS_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/GRS/")
 
#PDF/GRV/Makefile.am
 
# pkglib_LTLIBRARIES = libGRVSherpa.la
#GITTAG = PDF/GRV include ../../git.make libGRVSherpa_la
  set(PDF_GRV_SOURCES  grvphoton.f GRVph_Fortran_Interface.C)
#localinc_HEADERS = GRVph_Fortran_Interface.H
  list(TRANSFORM PDF_GRV_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/GRV/")
 
#PDF/Main/Makefile.am
 
# pkglib_LTLIBRARIES = libPDF.la
#GITTAG = PDF/Main include ../../git.make libPDF_la
  set(PDF_Main_SOURCES  PDF_Base.C ISR_Base.C Intact.C Structure_Function.C ISR_Handler.C Jet_Criterion.C Cluster_Definitions_Base.C NLOMC_Base.C Shower_Base.C Photon_PDF_Base.C)
#localinc_HEADERS = PDF_Base.H ISR_Base.H Intact.H Structure_Function.H\ ISR_Handler.H Jet_Criterion.H Cluster_Definitions_Base.H NLOMC_Base.H Shower_Base.H Photon_PDF_Base.H
  list(TRANSFORM PDF_Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/Main/")
 
#PDF/Makefile.am
 
# SUBDIRS = Main LHAPDF CT14 GRV SAL SASG GRS Electron NNPDF
 
#PDF/NNPDF/Makefile.am
 
# pkglib_LTLIBRARIES = libNNPDFSherpa.la
#GITTAG = PDF/NNPDF include ../../git.make libNNPDFSherpa_la
  set(PDF_NNPDF_SOURCES  NNPDFDriver.cc PDF_NN.C)
#localinc_HEADERS = NNPDFDriver.h NNPDF3archive = NNPDF31_nnlo_as_0118_mc.tar.bz2
#nobase_dist_pkgdata_DATA = $(NNPDF3archive) install-data-hook: tar xjf $(srcdir)/$(NNPDF3archive) -C $(DESTDIR)$(pkgdatadir) uninstall-hook: -rm -rf $(DESTDIR)/$(pkgdatadir)/NNPDF31_nnlo_as_0118_mc
  list(TRANSFORM PDF_NNPDF_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/NNPDF/")
 
#PDF/SAL/Makefile.am
 
# pkglib_LTLIBRARIES = libSALSherpa.la AM_CXXFLAGS = "-Wno-write-strings"
#GITTAG = PDF/SAL include ../../git.make libSALSherpa_la
  set(PDF_SAL_SOURCES  sal.cpp SALph_CPP_Interface.C)
#localinc_HEADERS = SALph_CPP_Interface.H sal.h math_ut.h
#nobase_dist_pkgdata_DATA = SALGrid/SAL.dat
  list(TRANSFORM PDF_SAL_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/SAL/")
 
#PDF/SASG/Makefile.am
 
# pkglib_LTLIBRARIES = libSASGSherpa.la
#GITTAG = PDF/SASG include ../../git.make libSASGSherpa_la
  set(PDF_SASG_SOURCES  sasgam2.f SASGph_Fortran_Interface.C)
#localinc_HEADERS = SASGph_Fortran_Interface.H
  list(TRANSFORM PDF_SASG_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../PDF/SASG/")
 
set(PDF_esources )
if (SHERPA_BUILD_STATIC)
add_library(PDF_static STATIC ${PDF_esources}
                             ${PDF_CT14_SOURCES} 
                             ${PDF_Electron_SOURCES} 
                             ${PDF_GRS_SOURCES} 
                             ${PDF_GRV_SOURCES} 
                             ${PDF_Main_SOURCES} 
                             ${PDF_NNPDF_SOURCES} 
                             ${PDF_SAL_SOURCES} 
                             ${PDF_SASG_SOURCES} 
)
target_include_directories(PDF_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PDF_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(PDF_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(PDF_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS PDF_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PDF_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PDF)
endif()
if (SHERPA_BUILD_SHARED)
add_library(PDF SHARED ${PDF_esources}
                             ${PDF_CT14_SOURCES} 
                             ${PDF_Electron_SOURCES} 
                             ${PDF_GRS_SOURCES} 
                             ${PDF_GRV_SOURCES} 
                             ${PDF_Main_SOURCES} 
                             ${PDF_NNPDF_SOURCES} 
                             ${PDF_SAL_SOURCES} 
                             ${PDF_SASG_SOURCES} 
)
target_include_directories(PDF PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PDF PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(PDF PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(PDF PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS PDF DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PDF        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PDF SOVERSION 0.0.0.${SHERPA_VERSION_MAJOR})
endif()
