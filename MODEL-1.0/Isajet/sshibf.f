CDECK  ID>, SSHIBF. 
      SUBROUTINE SSHIBF(IMODEL)
C-----------------------------------------------------------------------
C
C     This subroutine calculates the decay widths for decays of the
C     Higgs scalars present in the minimal SUSY model.
C
C     NOTE: Decays into sfermions are not yet incorporated.
C
C     Standard model parameters are hard wired in  SSMSSM. To get
C     the 1987-8 values corresponding to the Gunion et al. papers
C     (Intl. J. Mod. Phys. 2(4):1035; Nucl. Phys. B307:445) you must
C     change
C          ALFA3 = 0.12  --> 0.136
C          AMW   = 80.0  --> 81.3
C          AMZ   = 91.17 --> 92.7
C
C     2/9/91:
C     I've modified the program slightly.  The ALPHA3 = 0.12 value
C     above is the recent empirical value from LEP. Using the equation
C     from page 220 in Barger and Phillips yields ALPHA3 = 0.136.
C
C     10/1/92:
C     Now includes vertex corrections for triple Higgs couplings.
C     (See Kunszt and Zwirner, CERN-TH.6150/91 for all but HH-HC-HC
C     correction which is in our Higgs --> SUSY paper: Baer et al.
C     FSU-HEP-920630 or UH-511-749-92)
C
C     Bisset's HIGSBF
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          SM ident code definitions. These are standard ISAJET but
C          can be changed.
      INTEGER IDUP,IDDN,IDST,IDCH,IDBT,IDTP
      INTEGER IDNE,IDE,IDNM,IDMU,IDNT,IDTAU
      INTEGER IDGL,IDGM,IDW,IDZ,IDH
      PARAMETER (IDUP=1,IDDN=2,IDST=3,IDCH=4,IDBT=5,IDTP=6)
      PARAMETER (IDNE=11,IDE=12,IDNM=13,IDMU=14,IDNT=15,IDTAU=16)
      PARAMETER (IDGL=9,IDGM=10,IDW=80,IDZ=90,IDH=81)
C          SUSY ident code definitions. They are chosen to be similar
C          to those in versions < 6.50 but may be changed.
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      INTEGER IDTAUL,IDTAUR
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
      PARAMETER (IDTAUL=10016,IDTAUR=20016)
      INTEGER IMODEL
C
C          Decays into fermions
      CALL SSHFF(IMODEL)
C          Loop decays into photons and gluons
      CALL SSHGM
      CALL SSHGL
C          Decays into WW(*), ZZ(*)
      CALL SSHWW
C          Decays into neutralinos and charginos
      CALL SSHNN
      CALL SSHCC
C          Decays into other Higgs bosons
      CALL SSHHX
C          Decays to sfermions
      CALL SSHSF
C          Normalize branching ratios
C
      CALL SSNORM(ISHL)
      CALL SSNORM(ISHH)
      CALL SSNORM(ISHA)
      CALL SSNORM(ISHC)
C
      RETURN
      END
