c nopmax is the maximum number of PDF parameters one can use.
c noemax is the maximum number of PDF's in a list. 
c npfmax is the maximum number of functional parameters
c nofmax is the maximum number of functionals
c linemax is the maximum number of lines in pdf description text
c
c When changing the parameters make sure the 'make' recompiles all 
c affected *.f files
c
      integer noemax,nopmax,npfmax,nofmax,linemax
      parameter (noemax=1000,nopmax=30,npfmax=10,nofmax=10,linemax=20)
