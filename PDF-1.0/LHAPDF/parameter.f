      subroutine listPDF(imem,parm)
      implicit none
      include 'parmsetup.f'
      character*16 s1
      integer i,j,mem,imem,noe,nop,listN,listP,type
      real*8 parmL(0:noemax,nopmax),parm(nopmax)
      save listN,listP,parmL
*
      mem=imem
      if (mem.gt.listN) then
c         write(*,*) 'Maximum number of PDFs in list exceeded: ',
c     .              mem,' > ',listN
c         write(*,*) 'Returning most likely PDF'
         mem=0
      endif
      if (mem.lt.0) then
         write(*,*) 'Negative PDF member requested: ',mem
         write(*,*) 'Returning most likely PDF'
         mem=0
      endif
      do i=1,listP
         parm(i)=parmL(mem,i)
      enddo
      return
*
      entry nopPDF(nop)
      nop=listP
      return
*
      entry numberPDF(noe)
      noe=listN
      return
*
      entry InitListPDF
      type=-1
      read(1,*) s1,listN,listP
      if (index(s1,'list').eq.1) then
         type=1
         do i=0,listN
            read(1,*) (parmL(i,j),j=1,listP)
         enddo
      endif
      if (type.lt.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Unknown parameter list type ',s1
         stop
      endif
      return
*
      end
