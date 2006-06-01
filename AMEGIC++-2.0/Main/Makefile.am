pkglib_LTLIBRARIES = libAmegic.la

libAmegic_la_SOURCES = Helicity.C Topology.C Polarisation.C Pol_Info.C Point.C \
                       Process_Info.C Process_Base.C Single_Process.C \
                       Single_Process_MHV.C Process_Group.C \
                       All_Processes.C Full_Decay_Table.C All_Decays.C Amegic.C \
                       AMEGIC_SCT.C

pkginclude_HEADERS = Helicity.H Topology.H Polarisation.H Pol_Info.H Point.H \
                     Process_Info.H Process_Base.H Single_Process.H \
                     Single_Process_MHV.H Process_Group.H \
                     All_Processes.H Full_Decay_Table.H All_Decays.H Amegic.H \
                     AMEGIC_SCT.H

AM_CXXFLAGS   = @SHERPAFLAGS@

INCLUDES = @AMEGICINCS@ @EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ @BEAMINCS@ @PDFINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
