CDECK  ID>, SIGPL.  
        REAL FUNCTION SIGPL(P2,I,J,G1,G2,CHW)
C-----------------------------------------------------------------------
C          SIGPL: Chargino mass matrix correction
C     Taken from Damien M. Pierce, Jonathan A. Bagger, Konstantin T. Matchev,
C     Ren-jie Zhang, Nucl.Phys.B491:3-67,1997, hep-ph/9606211
C     Programmed by Tadas Krupovnickas
C          P2 = 4-momentum squared
C          CHW = Cos(theta_W) in DR bar scheme
C     Ordering: u=1,s=2,t=3,d=4,c=5,b=6,e=7,mu=8,tau=9,nue=10,num=11,nut=12
C     I and J are indexes in Pierce's base. To convert to ISAJET's basis
C     transform the matrix elements in the following way:
C     11->11, 12->-21, 21->-12, 22->22
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
C          MTQ                  = top mass at MSUSY
C          MBQ                  = bottom mass at MSUSY
C          MLQ                  = tau mass at MSUSY
C          FBMA                 = b-Yukawa at mA scale
C          VUQ                  = Hu vev at MSUSY
C          VDQ                  = Hd vev at MSUSY
      COMMON/SSPAR/AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,
     $VUQ,VDQ
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS,MTQ,MBQ,MLQ,FBMA,VUQ,VDQ
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C          Frozen couplings from RG equations:
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_h1^2     GSS(14) = M_h2^2     GSS(15) = M_er^2
C     GSS(16) = M_el^2     GSS(17) = M_dnr^2    GSS(18) = M_upr^2
C     GSS(19) = M_upl^2    GSS(20) = M_taur^2   GSS(21) = M_taul^2
C     GSS(22) = M_btr^2    GSS(23) = M_tpr^2    GSS(24) = M_tpl^2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = log(vdq)
C     GSS(31) = log(vuq)
C          Masses:
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C          Unification:
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUT
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      SAVE /SUGMG/
      COMMON/SSINF/XLAM
      DOUBLE PRECISION XLAM
      COMPLEX*16 SSB1,AC0PPW(4,2),AC0PPH(4,2,2)
      DOUBLE PRECISION TEMP,APPFF(2,12,2),APPCPZ(2,2)
      DOUBLE PRECISION APPCPG(2,2),APPCPH(2,2,4)
      REAL APPR(2,12),APPL(2,12),AP0PPW(4,2),APPPPZ(2,2)
     $,AP0PPH(4,2,2),APPPPH(2,2,4)
      REAL P2,G1,G2,CHW,COST,SINT,COSB,SINB,COSL,SINL,GG,GGP
      REAL THX,THY,VMAT(2,2),UMAT(2,2),BE,FB,FT,FMIX(12,2,2)
      INTEGER THIW1,THIW2
      INTEGER I,J,II,III
      COMPLEX IMAG,N(4,4)
      PARAMETER (IMAG=(0.,1.))
      REAL PI,SR2
      PI=4*ATAN(1.)
      SR2=SQRT(2.)
      COST=COS(THETAT)
      SINT=SIN(THETAT)
      COSB=COS(THETAB)
      SINB=SIN(THETAB)
      COSL=COS(THETAL)
      SINL=SIN(THETAL)
      GG=G2
      GGP=SQRT(3./5.)*G1
      BE=ATAN(VUQ/VDQ)
      FB=MBQ/VDQ
      FT=MTQ/VUQ
      THX=SIGN(1.,TAN(GAMMAL))
      THY=SIGN(1.,TAN(GAMMAR))
      IF (SIGN(1.,AMW1SS).EQ.1.) THEN
         THIW1=0
      ELSE
         THIW1=1
      END IF
      IF (SIGN(1.,AMW2SS).EQ.1.) THEN
         THIW2=0
      ELSE
         THIW2=1
      END IF
      VMAT(1,1)=-(-1)**THIW1*SIN(GAMMAR)
      VMAT(1,2)=(-1)**THIW1*COS(GAMMAR)
      VMAT(2,1)=-(-1)**THIW2*THY*COS(GAMMAR)
      VMAT(2,2)=-(-1)**THIW2*THY*SIN(GAMMAR)
      UMAT(1,1)=-SIN(GAMMAL)
      UMAT(1,2)=COS(GAMMAL)
      UMAT(2,1)=-THX*COS(GAMMAL)
      UMAT(2,2)=-THX*SIN(GAMMAL)
      DO II=1,4
        IF (SIGN(1.,AMZISS(II)).EQ.1.) THEN
          III=0
        ELSE
          III=1
        END IF
        N(II,1)=-IMAG**III*ZMIXSS(4,II)
        N(II,2)=-IMAG**III*ZMIXSS(3,II)
        N(II,3)=IMAG**III*ZMIXSS(2,II)
        N(II,4)=IMAG**III*ZMIXSS(1,II)
      ENDDO
      DO II=1,2
        DO III=1,12
          APPR(II,III)=0.
          APPL(II,III)=0.
        ENDDO
      ENDDO
      APPL(1,4)=GG
      APPL(1,5)=GG
      APPL(1,6)=GG
      APPL(1,7)=GG
      APPL(1,8)=GG
      APPL(1,9)=GG
      APPL(2,3)=-FT
      APPR(2,6)=-FT
      DO II=1,12
        FMIX(II,1,1)=1.
        FMIX(II,1,2)=0.
        FMIX(II,2,1)=0.
        FMIX(II,2,2)=1.
      ENDDO
      FMIX(3,1,1)=COST
      FMIX(3,1,2)=-SINT
      FMIX(3,2,1)=SINT
      FMIX(3,2,2)=COST
      FMIX(6,1,1)=COSB
      FMIX(6,1,2)=-SINB
      FMIX(6,2,1)=SINB
      FMIX(6,2,2)=COSB
      FMIX(9,1,1)=COSL
      FMIX(9,1,2)=-SINL
      FMIX(9,2,1)=SINL
      FMIX(9,2,2)=COSL
      DO II=1,2
        DO III=1,12
          APPFF(II,III,1)=FMIX(III,1,1)*APPL(II,III)
     $+FMIX(III,1,2)*APPR(II,III)
          APPFF(II,III,2)=FMIX(III,2,1)*APPL(II,III)
     $+FMIX(III,2,2)*APPR(II,III)
        ENDDO
      ENDDO
      DO II=1,4
        DO III=1,2
          AP0PPW(II,III)=0.
        ENDDO
      ENDDO
      AP0PPW(2,1)=-GG
      AP0PPW(4,2)=GG/SR2
      DO II=1,4
        DO III=1,2
          AC0PPW(II,III)=CONJG(N(II,1))*AP0PPW(1,III)
     $+CONJG(N(II,2))*AP0PPW(2,III)+CONJG(N(II,3))*AP0PPW(3,III)
     $+CONJG(N(II,4))*AP0PPW(4,III)
        ENDDO
      ENDDO
      DO II=1,2
        DO III=1,2
          APPPPZ(II,III)=0.
        ENDDO
      ENDDO
      APPPPZ(1,1)=GG*CHW
      APPPPZ(2,2)=GG*(2.*CHW**2-1.)/2./CHW
      DO II=1,2
        DO III=1,2
          APPCPZ(II,III)=VMAT(II,1)*APPPPZ(1,III)
     $+VMAT(II,2)*APPPPZ(2,III)
        ENDDO
      ENDDO
      DO II=1,2
        DO III=1,2
          APPCPG(II,III)=SQRT(4.*PI/137.036)*VMAT(III,II)
        ENDDO
      ENDDO
      DO II=1,4
        DO III=1,2
          AP0PPH(II,III,1)=0.
          AP0PPH(II,III,2)=0.
        ENDDO
      ENDDO
      AP0PPH(1,2,1)=GGP/SR2
      AP0PPH(2,2,1)=GG/SR2
      AP0PPH(3,1,1)=-GG
      DO II=1,4
        DO III=1,2
          AC0PPH(II,III,1)=COS(BE)*(CONJG(N(II,1))*AP0PPH(1,III,1)
     $+CONJG(N(II,2))*AP0PPH(2,III,1)+CONJG(N(II,3))*AP0PPH(3,III,1)
     $+CONJG(N(II,4))*AP0PPH(4,III,1))+SIN(BE)
     $*(CONJG(N(II,1))*AP0PPH(1,III,2)
     $+CONJG(N(II,2))*AP0PPH(2,III,2)+CONJG(N(II,3))*AP0PPH(3,III,2)
     $+CONJG(N(II,4))*AP0PPH(4,III,2))
          AC0PPH(II,III,2)=-SIN(BE)*(CONJG(N(II,1))*AP0PPH(1,III,1)
     $+CONJG(N(II,2))*AP0PPH(2,III,1)+CONJG(N(II,3))*AP0PPH(3,III,1)
     $+CONJG(N(II,4))*AP0PPH(4,III,1))+COS(BE)
     $*(CONJG(N(II,1))*AP0PPH(1,III,2)
     $+CONJG(N(II,2))*AP0PPH(2,III,2)+CONJG(N(II,3))*AP0PPH(3,III,2)
     $+CONJG(N(II,4))*AP0PPH(4,III,2))
        ENDDO
      ENDDO
      DO II=1,2
        DO III=1,2
          APPPPH(II,III,1)=0.
          APPPPH(II,III,2)=0.
          APPPPH(II,III,3)=0.
          APPPPH(II,III,4)=0.
        ENDDO
      ENDDO
      APPPPH(1,2,1)=GG/SR2
      APPPPH(2,1,2)=GG/SR2
      APPPPH(1,2,3)=GG/SR2
      APPPPH(2,1,4)=-GG/SR2
      DO II=1,2
        DO III=1,2
          APPCPH(II,III,1)=(UMAT(III,1)*APPPPH(II,1,1)
     $+UMAT(III,2)*APPPPH(II,2,1))*COS(ALFAH)
     $-(UMAT(III,1)*APPPPH(II,1,2)+UMAT(III,2)*APPPPH(II,2,2))
     $*SIN(ALFAH)
          APPCPH(II,III,2)=(UMAT(III,1)*APPPPH(II,1,1)
     $+UMAT(III,2)*APPPPH(II,2,1))*SIN(ALFAH)
     $+(UMAT(III,1)*APPPPH(II,1,2)+UMAT(III,2)*APPPPH(II,2,2))
     $*COS(ALFAH)
          APPCPH(II,III,3)=(UMAT(III,1)*APPPPH(II,1,3)
     $+UMAT(III,2)*APPPPH(II,2,3))*COS(BE)
     $+(UMAT(III,1)*APPPPH(II,1,4)+UMAT(III,2)*APPPPH(II,2,4))*SIN(BE)
          APPCPH(II,III,4)=-(UMAT(III,1)*APPPPH(II,1,3)
     $+UMAT(III,2)*APPPPH(II,2,3))*SIN(BE)
     $+(UMAT(III,1)*APPPPH(II,1,4)+UMAT(III,2)*APPPPH(II,2,4))*COS(BE)
        ENDDO
      ENDDO
      TEMP=DBLE((3.*(APPFF(I,1,1)*APPFF(J,1,1)
     $*SSB1(P2,AMUP,ABS(MSS(4)))
     $+APPFF(I,1,2)*APPFF(J,1,2)*SSB1(P2,AMUP,ABS(MSS(5)))
     $+APPFF(I,2,1)*APPFF(J,2,1)*SSB1(P2,AMST,ABS(MSS(8)))
     $+APPFF(I,2,2)*APPFF(J,2,2)*SSB1(P2,AMST,ABS(MSS(9)))
     $+APPFF(I,3,1)*APPFF(J,3,1)*SSB1(P2,AMTP,ABS(MSS(10)))
     $+APPFF(I,3,2)*APPFF(J,3,2)*SSB1(P2,AMTP,ABS(MSS(11)))
     $+APPFF(I,4,1)*APPFF(J,4,1)*SSB1(P2,AMDN,ABS(MSS(2)))
     $+APPFF(I,4,2)*APPFF(J,4,2)*SSB1(P2,AMDN,ABS(MSS(3)))
     $+APPFF(I,5,1)*APPFF(J,5,1)*SSB1(P2,AMCH,ABS(MSS(6)))
     $+APPFF(I,5,2)*APPFF(J,5,2)*SSB1(P2,AMCH,ABS(MSS(7)))
     $+APPFF(I,6,1)*APPFF(J,6,1)*SSB1(P2,AMBT,ABS(MSS(13)))
     $+APPFF(I,6,2)*APPFF(J,6,2)*SSB1(P2,AMBT,ABS(MSS(12))))
     $+APPFF(I,7,1)*APPFF(J,7,1)*SSB1(P2,AME,ABS(MSS(14)))
     $+APPFF(I,8,1)*APPFF(J,8,1)*SSB1(P2,AMMU,ABS(MSS(15)))
     $+APPFF(I,9,1)*APPFF(J,9,1)*SSB1(P2,AMTAU,ABS(MSS(16)))
     $+APPFF(I,10,1)*APPFF(J,10,1)*SSB1(P2,0.,ABS(MSS(17)))
     $+APPFF(I,10,2)*APPFF(J,10,2)*SSB1(P2,0.,ABS(MSS(18)))
     $+APPFF(I,11,1)*APPFF(J,11,1)*SSB1(P2,0.,ABS(MSS(19)))
     $+APPFF(I,11,2)*APPFF(J,11,2)*SSB1(P2,0.,ABS(MSS(20)))
     $+APPFF(I,12,1)*APPFF(J,12,1)*SSB1(P2,0.,ABS(MSS(21)))
     $+APPFF(I,12,2)*APPFF(J,12,2)*SSB1(P2,0.,ABS(MSS(22))))/2.
     $+CONJG(AC0PPW(1,I))*AC0PPW(1,J)*SSB1(P2,ABS(MSS(23)),AMW)
     $+CONJG(AC0PPW(2,I))*AC0PPW(2,J)*SSB1(P2,ABS(MSS(24)),AMW)
     $+CONJG(AC0PPW(3,I))*AC0PPW(3,J)*SSB1(P2,ABS(MSS(25)),AMW)
     $+CONJG(AC0PPW(4,I))*AC0PPW(4,J)*SSB1(P2,ABS(MSS(26)),AMW)
     $+APPCPZ(I,1)*APPCPZ(J,1)*SSB1(P2,ABS(MSS(27)),AMZ)
     $+APPCPZ(I,2)*APPCPZ(J,2)*SSB1(P2,ABS(MSS(28)),AMZ)
     $+APPCPG(I,1)*APPCPG(J,1)*SSB1(P2,ABS(MSS(27)),0.)
     $+APPCPG(I,2)*APPCPG(J,2)*SSB1(P2,ABS(MSS(28)),0.)
     $+(CONJG(AC0PPH(1,I,1))*AC0PPH(1,J,1)
     $*SSB1(P2,ABS(MSS(23)),AMW)
     $+CONJG(AC0PPH(1,I,2))*AC0PPH(1,J,2)
     $*SSB1(P2,ABS(MSS(23)),ABS(MSS(32)))
     $+CONJG(AC0PPH(2,I,1))*AC0PPH(2,J,1)*SSB1(P2,ABS(MSS(24)),AMW)
     $+CONJG(AC0PPH(2,I,2))*AC0PPH(2,J,2)
     $*SSB1(P2,ABS(MSS(24)),ABS(MSS(32)))
     $+CONJG(AC0PPH(3,I,1))*AC0PPH(3,J,1)*SSB1(P2,ABS(MSS(25)),AMW)
     $+CONJG(AC0PPH(3,I,2))*AC0PPH(3,J,2)
     $*SSB1(P2,ABS(MSS(25)),ABS(MSS(32)))
     $+CONJG(AC0PPH(4,I,1))*AC0PPH(4,J,1)*SSB1(P2,ABS(MSS(26)),AMW)
     $+CONJG(AC0PPH(4,I,2))*AC0PPH(4,J,2)
     $*SSB1(P2,ABS(MSS(26)),ABS(MSS(32))))/2.
     $+(APPCPH(I,1,1)*APPCPH(J,1,1)
     $*SSB1(P2,ABS(MSS(27)),ABS(MSS(30)))
     $+APPCPH(I,1,2)*APPCPH(J,1,2)
     $*SSB1(P2,ABS(MSS(27)),ABS(MSS(29)))
     $+APPCPH(I,1,3)*APPCPH(J,1,3)*SSB1(P2,ABS(MSS(27)),AMZ)
     $+APPCPH(I,1,4)*APPCPH(J,1,4)
     $*SSB1(P2,ABS(MSS(27)),ABS(MSS(31)))
     $+APPCPH(I,2,1)*APPCPH(J,2,1)
     $*SSB1(P2,ABS(MSS(28)),ABS(MSS(30)))
     $+APPCPH(I,2,2)*APPCPH(J,2,2)
     $*SSB1(P2,ABS(MSS(28)),ABS(MSS(29)))
     $+APPCPH(I,2,3)*APPCPH(J,2,3)*SSB1(P2,ABS(MSS(28)),AMZ)
     $+APPCPH(I,2,4)*APPCPH(J,2,4)
     $*SSB1(P2,ABS(MSS(28)),ABS(MSS(31))))/2.
     $)/16./PI**2
      SIGPL=TEMP
      RETURN
      END
