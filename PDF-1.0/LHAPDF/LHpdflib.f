      subroutine InitPDFset(name)
      implicit none
      character name*(*)
      character*64 string
      character*16 s1,s2
      integer id,token,Ctoken
*
      open(unit=1,file=name,status='old')
      read(1,*) s1,s2
      if ((index(s2,'1.0').ne.1).and.(index(s2,'1.1').ne.1)) then
         write(*,*) 
     .        'Version ',s2,' not supported by this version of LHPDF'
         stop
      else  
         write(*,*) '**************************************************'
         write(*,*) '* LHAPDF Version 1.0 release                     *'
         write(*,*) '**************************************************'
         write(*,*)
      endif
      id=Ctoken()
 1    read(1,*) string
      id=token(string)
      if (id.eq.0) then
         write(*,*) 'File description error:'
         write(*,*) 'Command not understood: ',string
         stop
      endif
      if (id.eq.1) call descriptionPDF(id)
      if (id.eq.2) call initEvolve
      if (id.eq.3) call initAlphasPDF
      if (id.eq.4) call initInputPDF
      if (id.eq.5) call initListPDF
      if (id.ne.6) goto 1
      close(1)
      call InitEvolveCode
*
      return
      end
*     
      integer function token(s)
      implicit none
      character*16 s
      integer not,i,Ctoken
      parameter(not=6)
      character*16 t(not)/'Description:','Evolution:','Alphas:',
     .                    'Parametrization:','Parameterlist:',
     .                    'End:'/
      integer count(not)
      save count
*
      token=0
      do i=1,not
         if (s.eq.t(i)) token=i
      enddo
      if (token.ne.0) then
         count(token)=count(token)+1
         if (count(token).eq.2) then
            write(*,*) 'File description error:'
            write(*,*) 'Second definition of entry: ',s
            stop
         endif
      endif
      return
*
      entry Ctoken()
      do i=1,not
         count(i)=0
      enddo
      Ctoken=0
      return
*     
      end
