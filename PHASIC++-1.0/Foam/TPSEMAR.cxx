#include"TPSEMAR.h"

#ifdef ROOT_DEF
ClassImp(TPSEMAR)
#endif
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//       PseudoRandom number generator used by Foam 2.x                    //
// Fortran version rewritten by S. Jadach, Nov 1997.                       //
// C++ translation by M. Slusarczyk, Feb 2000, and S. Jadach, Oct 2000	   // 
//                                                                         //
// Universal random number generator proposed by MARSAGLIA and ZAMAN       //
// in report FSU-SCRI-87-50                                                //
//        modified by F. James, 1988 and 1989, to generate a vector        //
//        of pseudorandom numbers rvec of length lenv, and to put in       //
//        the COMMON block everything needed to specify currrent state,    //
//        and to add input and output entry points rmarin, rmarut.         //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   //
// ++  CALLing sequences for TPSEMAR:                                 ++   //
// ++      TPSEMAR()          INITIALIZES the generator               ++   //
// ++      TPSEMAR(i1,n1,n2)  INITIALIZES the generator               ++   //
// ++                  32-bit integer i1, and number counts n1,n2     ++   //
// ++                  (for initializing, set n1=n2=0, but to restart ++   //
// ++                    a previously generated sequence, use values  ++   //
// ++                    output by rmarut)                            ++   //
// ++      Initialize(i1,n1,n2) used by the two above constructors    ++   //
// ++      MakeVec(lenv,RVec) GENERATES random vector RVec            ++   //
// ++      Out(i1,n1,n2)   outputs the value of the original          ++   //
// ++                  seed and the two number counts, to be used     ++   //
// ++                  for restarting by initializing to i1 and       ++   //
// ++                  skipping n2*100000000+n1 numbers.              ++   //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   //
//                                                                         //
//         Initializing routine for ranmar, may be called before           //
//         generating pseudorandom numbers with ranmar. the input          //
//         values should be in the ranges:  0<=ijklin<=900 000 000         //
//                                          0<=ntotin<=999 999 999         //
//                                          0<=ntot2n<<999 999 999!        //
// to get the standard values in MARSAGLIA's paper, ijklin=54217137        //
//                                            ntotin,ntot2n=0              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
TPSEMAR::TPSEMAR(void){
  //   Constructor with internal seed
  m_modcns = 1000000000;
  long ijkl_new  = 54217137;
  long ntot_new  = 0;
  long ntot2_new = 0;
  Initialize(ijkl_new, ntot_new,ntot2_new);
}// TPSEMAR
/////////////////////////////////////////////////////////////////////////////
TPSEMAR::TPSEMAR(long ijkl_new, long ntot_new, long ntot2_new){
  //   Constructor with external seed
  Initialize(ijkl_new, ntot_new,ntot2_new);
}
/////////////////////////////////////////////////////////////////////////////
TPSEMAR::~TPSEMAR(){
  //   Destructor
}
/////////////////////////////////////////////////////////////////////////////
void TPSEMAR::Initialize(long ijkl_new, long ntot_new, long ntot2_new){
  // RE-initialize random number generator
  // prividing seed numbers
  double t,uni,s;
  long int m,i24,jj,idum,loop2,now;
  long int i,j,ij,k,l,ii,kl;
  
  m_ijkl = ijkl_new;
  m_ntot = max(ntot_new,0);
  m_ntot2= max(ntot2_new,0);
  
  ij = m_ijkl/30082;
  kl = m_ijkl - 30082*ij;
  i = MOD(ij/177, 177) + 2;
  j = MOD(ij, 177)     + 2;
  k = MOD(kl/169, 178) + 1;
  l = MOD(kl, 169);
  
  for( ii= 1; ii<=97 ; ii++)
    {
      s = 0.0;
      t = 0.5;
      for(jj= 1; jj<= 24 ;jj++)
	{
	  m = MOD(MOD(i*j,179)*k, 179);
	  i = j;
	  j = k;
	  k = m;
	  l = MOD(53*l+1, 169);
	  if (MOD(l*m,64) >=  32)  s = s+t;
	  t = 0.5*t;
	}
      m_U[ii] = s;
    }
  
  m_twom24 = 1.0;
  
  for(i24 = 1 ; i24 <= 24 ; i24++)
    m_twom24 = 0.5*m_twom24;
  
  m_C  =   362436.*m_twom24;
  m_CD =  7654321.*m_twom24;
  m_CM = 16777213.*m_twom24;
  m_i97 = 97;
  m_j97 = 33;
  
  // complete initialization by skipping
  // (ntot2*m_modcns + ntot) random numbers  
  for(loop2= 1; loop2<= m_ntot2+1 ; loop2++)
    {
      now = m_modcns;
      if (loop2  ==  m_ntot2+1)  now=m_ntot;
      if (now  >  0)
	{
	  for(idum = 1; idum <= m_ntot ; idum++)
	    {
	      uni = m_U[m_i97]-m_U[m_j97];
	      if (uni  <  0.)  uni=uni+1.0;
	      m_U[m_i97] = uni;
	      m_i97 = m_i97-1;
	      if (m_i97  ==  0)  m_i97=97;
	      m_j97 = m_j97-1;
	      if (m_j97  ==  0)  m_j97=97;
	      m_C = m_C - m_CD;
	      if (m_C  <  0.)   m_C=m_C+m_CM;
	    }
	}
    }
}// Initialize
/////////////////////////////////////////////////////////////////////////////
void TPSEMAR::MakeVec(int lenv, double *RVec){
  //  Generate lenv random numbers
  double   zuni,uni;
  int ivec;
  // Normal entry to generate lenv random numbers
  for( ivec= 0; ivec < lenv; ivec++)
    {
      uni = m_U[m_i97]-m_U[m_j97];
      if (uni  <  0.)  uni=uni+1.0;
      m_U[m_i97] = uni;
      m_i97 = m_i97-1;
      if (m_i97  ==  0)  m_i97=97;
      m_j97 = m_j97-1;
      if (m_j97  ==  0)  m_j97=97;
      m_C = m_C - m_CD;
      if (m_C  <  0.)   m_C=m_C+m_CM;
      uni = uni-m_C;
      if (uni  <  0.) uni=uni+1.0;
      *(RVec +ivec) = uni;
      // Replace exact zeros by uniform distr. *2**-24
      if (uni  ==  0.)
	{
	  zuni = m_twom24*m_U[2];
	  // an exact zero here is very unlikely, but let's be safe.
	  if (zuni  ==  0.) zuni= m_twom24*m_twom24;
	  *(RVec +ivec) = zuni;
	}
    }
  m_ntot  = m_ntot + lenv;
  if (m_ntot  >=  m_modcns)
    {
      m_ntot2  =  m_ntot2 + 1;
      m_ntot   =  m_ntot - m_modcns;
    }
}//MakeVec
/////////////////////////////////////////////////////////////////////////////
void TPSEMAR::Out(long int &ijkl_out,long int &ntot_out,long int &ntot2_out){
  // Get position of the rn generator for later restart
  ijkl_out  = m_ijkl;
  ntot_out  = m_ntot;
  ntot2_out = m_ntot2;
}
/////////////////////////////////////////////////////////////////////////////
//          End of class  TPSEMAR                                          //
/////////////////////////////////////////////////////////////////////////////
