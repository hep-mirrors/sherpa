#include "Random.H"
#include "Message.H"
#include "MathTools.H"
#include "MyStrStream.H"
#include <iostream>

using namespace ATOOLS;
using namespace std;

#define MAXLOGFILES 10

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

ATOOLS::Random ATOOLS::ran(1234);

Random::Random(long nid): 
  p_outstream(NULL) 
{ 
  SetSeed(nid); 
  Exception_Handler::AddTerminatorObject(this);
  SaveStatus();
}

Random::~Random() 
{ 
  if (p_outstream!=NULL) {
    p_outstream->close();
    delete p_outstream;
  }
} 

static long idum2=123456789;
static long sidum2=123456789;
static long iy=0;
static long siy=0;
static long iv[NTAB];
static long siv[NTAB];

double Random::Ran2(long *idum)
{
  int   j;
  long  k;
  double temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Random::Ran1(long *idum)
{
  int j;
  long k;
  double temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

void Random::InitRan3(long *idum)
{
  long mj,mk;
  int i,ii,k;
  
  mj=MSEED-(*idum < 0 ? -*idum : *idum);
  mj %= MBIG;
  m_ma[55]=mj;
  mk=1;
  for (i=1;i<=54;i++) {
    ii=(21*i) % 55;
    m_ma[ii]=mk;
    mk=mj-mk;
    if (mk < MZ) mk += MBIG;
    mj=m_ma[ii];
  }
  for (k=1;k<=4;k++)
    for (i=1;i<=55;i++) {
      m_ma[i] -= m_ma[1+(i+30) % 55];
      if (m_ma[i] < MZ) m_ma[i] += MBIG;
    }
  m_inext=0;
  m_inextp=31;
}

double Random::Ran3()
{
  if (++m_inext == 56) m_inext=1;
  if (++m_inextp == 56) m_inextp=1;
  long mj=m_ma[m_inext]-m_ma[m_inextp];
  if (mj < MZ) mj += MBIG;
  m_ma[m_inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */


int Random::WriteOutStatus(const char * filename){
  // write out every Statusregister of Random Number generator

  //  sprintf(m_outname,"%s%i.dat",filename,m_written); 
  if ((p_outstream!=0) && (std::strcmp(filename,m_outname)!=0)) {
    p_outstream->close();
    p_outstream = 0;
  }
  if (p_outstream == 0){
    msg_Tracking()<<"Random::WriteOutStatus : Saving Random Number Generator Status to "<<filename<<endl;
    long int count=0;
    std::ifstream *myinstream = new std::ifstream(filename,std::ios::in);
    if (myinstream->good()) {
      std::string buffer;
      while (!myinstream->eof()) {
	(*myinstream)>>count;
	getline(*myinstream,buffer);
      }
      ++count;
      myinstream->close();
      delete myinstream;
    }
#ifdef _IOS_BAD
    p_outstream = new std::fstream(filename,ios::app | ios::out);
#else
    p_outstream = new std::fstream(filename,std::ios_base::app | std::ios_base::out);
#endif
    std::strcpy(m_outname,filename);
    m_written=count;
  } 
  (*p_outstream)<<m_written<<"\t"<<m_id<<"\t"<<m_inext<<"\t"<<m_inextp<<"\t";
  for (int i=0;i<56;++i) (*p_outstream)<<m_ma[i]<<"\t";
  (*p_outstream)<<iy<<"\t"<<idum2<<"\t";
  for (int i=0;i<NTAB;++i) (*p_outstream)<<iv[i]<<"\t";
  (*p_outstream)<<endl;
  return m_written++;
}

void Random::ReadInStatus(const char * filename, long int index){
  // read in every Statusregister of Random Number generator
  msg_Info()<<"Random::ReadInStatus from "<<filename<<" index "<<index<<endl;
  std::ifstream myinstream(filename);
  long int count;
  if (myinstream.good()) {
    (myinstream)>>count;
    std::string buffer;
    while (count!=index && !myinstream.eof()) {
      getline(myinstream,buffer);
      (myinstream)>>count;    
    }
    if (count==index) {
      (myinstream)>>m_id; (myinstream)>>m_inext; (myinstream)>>m_inextp;
      for (int i=0;i<56;++i) (myinstream)>>m_ma[i];    
      (myinstream)>>iy>>idum2;
      for (int i=0;i<NTAB;++i) (myinstream)>>iv[i];
    } 
    else msg.Error()<<"ERROR in Random::ReadInStatus : index="<<index<<" not found in "<<filename<<endl;
    myinstream.close();
  } 
  else msg.Error()<<"ERROR in Random::ReadInStatus : "<<filename<<" not found!!"<<endl;
}

double Random::GetNZ() 
{
  double ran1;
  do ran1=Get(); while (ran1==0.); 
  return ran1;
}

void Random::SetSeed(long int nid) 
{
  m_id = nid<0 ? nid : -nid;
  InitRan3(&m_id);
  m_written=0;    
  p_outstream=0;
  std::strcpy(m_outname,"");
}

void Random::SaveStatus()
{
  m_sid=m_id; 
  m_sinext=m_inext; 
  m_sinextp=m_inextp;
  for (int i=0;i<56;++i) m_sma[i]=m_ma[i];    
  siy=iy;
  sidum2=idum2;
  for (int i=0;i<NTAB;++i) siv[i]=iv[i];
}

void Random::RestoreStatus()
{
  m_id=m_sid; 
  m_inext=m_sinext; 
  m_inextp=m_sinextp;
  for (int i=0;i<56;++i) m_ma[i]=m_sma[i];    
  iy=siy;
  idum2=sidum2;
  for (int i=0;i<NTAB;++i) iv[i]=siv[i];
}

void Random::PrepareTerminate()
{
  if (Exception_Handler::LastException()==NULL && 
      Exception_Handler::LastSignal()==0) return;
  if (Exception_Handler::LastException()!=NULL) {
    if (Exception_Handler::LastException()->Type()==ex::normal_exit) return;
  }
  if (Exception_Handler::LastSignal()==2) return;
  std::string name="debug_info_";
  unsigned int i=0;
  do { 
    std::ifstream testfile((name+ATOOLS::ToString(++i)+
			    std::string(".random")).c_str());
    if (!testfile.is_open()) break;
  } while (i<MAXLOGFILES);
  RestoreStatus();
  WriteOutStatus((name+ATOOLS::ToString(i)+
		  std::string(".random")).c_str());
}
