      subroutine descriptionPDF(id)
      implicit none
      include 'parmsetup.f'
      integer id,token,l,nline
      character*64 string
      character*64 desc(linemax)
      save nline,desc
*
      l=0
      write(*,*) '>>>>>> PDF description: <<<<<<'
 1    read(1,*) string
      id=token(string)
      if (id.eq.0) then
         write(*,*) string
         l=l+1
         if (l.gt.linemax) then
            write(*,*) 'Too many lines in PDF description to store.'
            write(*,*) 'Increase linemax variable in parmsetup.f'
            write(*,*) 'Ignoring additional description lines'
            l=linemax
         else
            desc(l)=string
            goto 1
         endif
      endif
      nline=l
      write(*,*)'>>>>>>                   <<<<<<'
      write(*,*)
      return
*
      entry GetDesc
      do l=1,nline
         write(*,*) desc(l)
      enddo
      return
*     
      end
