C******************************************************************
C
C     Wrappers around CTEQ6 routines
C
C******************************************************************

C******************************************************************

      SUBROUTINE CTQ6INITSET(set)
      integer set
C      print *,set
      call SetCtq6(set)
      end

C******************************************************************

      FUNCTION CTQ6EVOLVE(f,x,Q)
      implicit none
      integer f
      real*8 x,Q
      double precision c1, CTQ6EVOLVE, Ctq6Pdf
C      print *,f
C      print *,'x = ',x
C      print *,'Q = ',Q
      CTQ6EVOLVE  = Ctq6Pdf (f,x,Q)
C      print *,'CTQ6EVOLVE = ',CTQ6EVOLVE
C      print *,'CTQ6EVOLVE = ',c1
      
      end

C******************************************************************
