C******************************************************************
C
C     Wrappers around LHAPDF routines
C
C******************************************************************

C******************************************************************

      SUBROUTINE LHAPDFINITSET(set)
      character*64 set

      print *,set
      call initpdfset(set)

      end

C******************************************************************

      SUBROUTINE LHAPDFINIT(member)
      integer member
      call initpdf(member)
      end

C******************************************************************

      SUBROUTINE LHAPDFEVOLVE(x,Q,f)
      real*8 x,Q
      real*8 f(0:12)
      call evolvepdf(x,Q,f)
      end

C******************************************************************

      FUNCTION LHAPDFALPHAS(Q)
      real*8 Q
      ALPHAS = alphasPDF(Q)
      end

C******************************************************************
