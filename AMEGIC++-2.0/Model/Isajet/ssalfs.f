CDECK  ID>, SSALFS.
      DOUBLE PRECISION FUNCTION SSALFS(Q2)
C-----------------------------------------------------------------------
C     Strong coupling formula from page 201 of Barger and Phillips:
C     (using ALQCD4 for 4 flavor Lambda)
C
C     Bisset's STRCPLH
C-----------------------------------------------------------------------
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
C
      DOUBLE PRECISION Q2,AS,TH5,TH6,PI
      DATA PI/3.14159265D0/
C
      TH5=4*AMBT**2
      TH6=4*AMTP**2
      IF (Q2.LE.TH5)THEN
        AS=12*PI/(25*LOG(Q2/ALQCD4**2))
      ELSE IF(Q2.GT.TH5.AND.Q2.LE.TH6) THEN
        AS=25*LOG(Q2/ALQCD4**2)-2*LOG(Q2/TH5)
        AS=12*PI/AS
      ELSEIF(Q2.GT.TH6)THEN
        AS=25*LOG(Q2/ALQCD4**2)
        AS=AS-2*(LOG(Q2/TH5)+LOG(Q2/TH6))
        AS=12*PI/AS
      ENDIF
      SSALFS=AS
      RETURN
      END
