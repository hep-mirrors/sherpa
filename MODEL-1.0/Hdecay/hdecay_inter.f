C =======================================================================
C
      SUBROUTINE HDECAYINTER(hmass,iflags,acouplings,
     .                       abosons,ayukawas,ackms)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DIMENSION iflags(5)
      DIMENSION acouplings(4),abosons(5),ayukawas(10),ackms(4)


      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/FLAGS_HDEC/INDIDEC
      PARAMETER(K=6,NI=87,NSA=85,NSB=86,NLA=88,NLB=89,NHA=90,NHB=91,
     .          NHC=92,NAA=93,NAB=94,NCA=95,NCB=96,NRA=97,NRB=98,
     .          NSUSYL=81,NSUSYA=82,NSUSYH=83,NSUSYC=84,NPAR=80,
     .          NSUSYLA=79,NSUSYLB=78,NSUSYLC=77,NSUSYLD=76,NSUSYLE=75,
     .          NSUSYLF=59,NSUSYHF=58,
     .          NSUSYHA=74,NSUSYHB=73,NSUSYHC=72,NSUSYHD=71,NSUSYHE=70,
     .          NSUSYAA=69,NSUSYAB=68,NSUSYAC=67,NSUSYAD=66,NSUSYAE=65,
     .          NSUSYCA=64,NSUSYCB=63,NSUSYCC=62,NSUSYCD=61,NSUSYCE=60)
      DIMENSION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .          GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .          HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4)
      DIMENSION GLTT(2,2),GLBB(2,2),GHTT(2,2),GHBB(2,2),GCTB(2,2),
     .          GLEE(2,2),GHEE(2,2),GCEN(2,2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/SFER1ST_HDEC/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,
     .              ABRZGA,ABRZ,AWDTH
      COMMON/WIDTHHL_HDEC/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,
     .               HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,
     .               HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,
     .               HHBRGA,HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,
     .               HHBRHW,HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,
     .               HCBRW,HCBRA,HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER_HDEC/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,
     .              BHLSQDL,BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS_HDEC/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD

      nma      = 1
      tgbet    = 1.d0
      ihiggs   = 0
      imodel   = 0
      ionsh    = 0

      nnlo     = iflags(1)  
      iowz     = iflags(2)
      ipole    = iflags(3)
      nfgg     = iflags(4)

      amabeg   = hmass
      amaend   = hmass
      alph     = acouplings(1)
      GF       = acouplings(2)
      alsmz    = acouplings(3)
      amW      = abosons(1)
      amZ      = abosons(2)
      GamW     = abosons(3)
      GamZ     = abosons(4)
      ams      = ayukawas(3)
      amc      = ayukawas(4)
      amb      = ayukawas(5)
      amt      = ayukawas(6)
      ammuon   = ayukawas(8)
      amtau    = ayukawas(9)
      VUS      = ackms(1)
      VCB      = ackms(2)
      VUB      = ackms(3)
      RVUB     = ackms(3)/ackms(2)
      if (VCB.lt.1.e-6.or.vcb.gt.-1.e-6) RVUB = 0.d0

      CALL HEAD_HDEC(TGBET,AMABEG)

      PI      = 4*DATAN(1D0)

      ALPH    = 1.D0/ALPH
      AMSB    = AMS

      AMC0    = AMC
      AMB0    = AMB
      AMT0    = AMT
      ACC     = 1.D-8
      NLOOP   = 2
      XLAMBDA = XITLA_HDEC(NLOOP,ALSMZ,ACC)
      N0=5
      CALL ALSINI_HDEC(ACC)

C--INITIALIZE COEFFICIENTS FOR POLYLOGARITHMS
      NBER = 18
      CALL BERNINI_HDEC(NBER)

C--CHECK NFGG
      IF(NFGG.GT.5.OR.NFGG.LT.3)THEN
       WRITE(6,*)'NF-GG NOT VALID. TAKING THE DEFAULT NF-GG = 3....'
       NFGG = 3
      ENDIF

100   FORMAT(10X,G30.20)
101   FORMAT(10X,I30)

C--WRITE THE INPUT PARAMTERS TO A DATA-FILE

      WRITE(NPAR,8)'HIGGS    = ',IHIGGS
      WRITE(NPAR,8)'MODEL    = ',IMODEL
      WRITE(NPAR,9)'TGBET    = ',TGBET
      WRITE(NPAR,9)'MABEG    = ',AMABEG
      WRITE(NPAR,9)'MAEND    = ',AMAEND
      WRITE(NPAR,7)'NMA      = ',NMA
      WRITE(NPAR,9)'ALS(MZ)  = ',ALSMZ
      WRITE(NPAR,9)'MSBAR(1) = ',AMS
      WRITE(NPAR,9)'MC       = ',AMC
      WRITE(NPAR,9)'MB       = ',AMB
      WRITE(NPAR,9)'MT       = ',AMT
      WRITE(NPAR,9)'MTAU     = ',AMTAU
      WRITE(NPAR,9)'MMUON    = ',AMMUON
      WRITE(NPAR,9)'ALPH     = ',1.D0/ALPH
      WRITE(NPAR,9)'GF       = ',GF
      WRITE(NPAR,9)'GAMW     = ',GAMW
      WRITE(NPAR,9)'GAMZ     = ',GAMZ
      WRITE(NPAR,9)'MZ       = ',AMZ
      WRITE(NPAR,9)'MW       = ',AMW
      WRITE(NPAR,9)'VUS      = ',VUS
      WRITE(NPAR,9)'VCB      = ',VCB
      WRITE(NPAR,9)'VUB/VCB  = ',RVUB
      WRITE(NPAR,9)'MU       = ',AMU
      WRITE(NPAR,9)'M2       = ',AM2
      WRITE(NPAR,9)'MSL1      = ',AMEL1
      WRITE(NPAR,9)'MER1      = ',AMER1
      WRITE(NPAR,9)'MQL1      = ',AMQL1
      WRITE(NPAR,9)'MUR1      = ',AMUR1
      WRITE(NPAR,9)'MDR1      = ',AMDR1
      WRITE(NPAR,9)'MEL      = ',AMEL
      WRITE(NPAR,9)'MER      = ',AMER
      WRITE(NPAR,9)'MSQ      = ',AMSQ
      WRITE(NPAR,9)'MUR      = ',AMUR
      WRITE(NPAR,9)'MDR      = ',AMDR
      WRITE(NPAR,9)'AL       = ',AL
      WRITE(NPAR,9)'AU       = ',AU
      WRITE(NPAR,9)'AD       = ',AD
      WRITE(NPAR,8)'NNLO (M) = ',NNLO
      WRITE(NPAR,8)'ON-SHELL = ',IONSH
      WRITE(NPAR,8)'ON-SH-WZ = ',IONWZ
      WRITE(NPAR,8)'OFF-SUSY = ',IOFSUSY
      WRITE(NPAR,8)'IPOLE    = ',IPOLE 
      WRITE(NPAR,8)'NF-GG    = ',NFGG
      WRITE(NPAR,9)'LAMBDA_5 = ',XLAMBDA

      CLOSE(NPAR)

7     FORMAT(A11,I7)
8     FORMAT(A11,I4)
9     FORMAT(A11,G15.6)

      CLOSE(NI)


      AMAR = AMABEG
      AMSM = AMAR
      AMA  = AMAR


      CALL HDEC(TGBET)
      CALL WRITE_HDEC(TGBET)

      CALL CLOSE_HDEC

      END

C =======================================================================
C =======================================================================
C =======================================================================


      SUBROUTINE HDECAYSM(BRSMFF,BRSMVV,WIDTH)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

      DIMENSION BRSMFF(10),BRSMVV(6)
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH

      BRSMFF(1) = 0.D0
      BRSMFF(2) = 0.D0
      BRSMFF(3) = smbrs
      BRSMFF(4) = smbrc
      BRSMFF(5) = smbrb
      BRSMFF(6) = smbrt
      BRSMFF(7) = 0.D0
      BRSMFF(8) = smbrm
      BRSMFF(9) = smbrl

      BRSMVV(1) = smbrg
      BRSMVV(2) = smbrga
      BRSMVV(3) = smbrw
      BRSMVV(4) = smbrz
      BRSMVV(5) = smbrzga

      WIDTH     = smwdth

      END
