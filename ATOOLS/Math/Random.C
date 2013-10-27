#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Marsaglia.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_MPI.H"
#include <iostream>
#include <cstring>

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::External_RNG
#define PARAMETER_TYPE ATOOLS::RNG_Key
#include "ATOOLS/Org/Getter_Function.C"

using namespace std;
using namespace ATOOLS;

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

ATOOLS::Random *ATOOLS::ran(NULL);


ATOOLS::Random::Random(long nid): 
  m_nsinceinit(0), m_increment(0), p_external(NULL)
{
  ATOOLS::exh->AddTerminatorObject(this);
  SetSeed(nid); 
  SaveStatus();
  p_ran4[0] = new Marsaglia();
  p_ran4[1] = new Marsaglia();
}


ATOOLS::Random::~Random() 
{ 
  delete p_ran4[0];
  delete p_ran4[1];
  if (p_external) delete p_external;
} 

static long idum2=123456789;
static long sidum2=123456789;
static long iy=0;
static long siy=0;
static long iv[NTAB];
static long siv[NTAB];

double ATOOLS::Random::Ran2(long *idum)
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


int ATOOLS::Random::WriteOutStatus(const char * filename)
{
  if (p_external!=NULL) return 1;
  // if Ran4 is the active Random Generator, then use its method
  if (activeGenerator==4) {return WriteOutStatus4(filename);}
  // write out every Statusregister of Random Number generator

  // remove old random file
  if (FileExists(filename)) remove(filename);
  std::ofstream outstream(filename);
  outstream<<0<<"\t"<<m_id<<"\t";
  outstream<<iy<<"\t"<<idum2<<"\t";
  for (int i=0;i<NTAB;++i) outstream<<iv[i]<<"\t";
  outstream<<endl;
  return 1;
}

int ATOOLS::Random::WriteOutSavedStatus(const char * filename)
{
  if (p_external!=NULL) return 1;
  // if Ran4 is the active Random Generator, then use its method
  if (activeGenerator==4) {return WriteOutSavedStatus4(filename);}
  // write out every Statusregister of Random Number generator

  // remove old random file
  if (FileExists(filename)) remove(filename);
  std::ofstream outstream(filename);
  outstream<<0<<"\t"<<m_sid<<"\t";
  outstream<<siy<<"\t"<<sidum2<<"\t";
  for (int i=0;i<NTAB;++i) outstream<<siv[i]<<"\t";
  outstream<<endl;
  return 1;
}

int ATOOLS::Random::WriteOutStatus
(std::ostream &outstream,const size_t &idx)
{
  if (activeGenerator==4) return WriteOutStatus4(outstream,idx);
  outstream<<idx<<"\t"<<m_id<<"\t";
  outstream<<iy<<"\t"<<idum2<<"\t";
  for (int i=0;i<NTAB;++i) outstream<<iv[i]<<"\t";
  outstream<<"\n";
  return idx+1;
}

bool ATOOLS::Random::ReadInStatus(const std::string &path) 
{
  ReadInStatus((path+"random.dat").c_str());
  return true;
}

void ATOOLS::Random::ReadInStatus(const char * filename)
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_size()>1) return;
#endif
  if (p_external!=NULL) return;
  // check what type of data is in target file
  ifstream file(filename);
  // Check if the first 20 bytes in file are a known identifier and
  // set the activeGenerator variable accordingly
  if (FileExists(std::string(filename)+".msg")) activeGenerator = 4;
  else activeGenerator = 2;
  file.close();
  
  // use readin method for the Generator identified Generator
  if (activeGenerator==4) { ReadInStatus4(filename);} 
  else {  
    // read in every Statusregister of Random Number generator
    msg_Info()<<METHOD<<"(): Reading status from '"
	      <<filename<<"'."<<endl;
    std::ifstream myinstream(filename);
    long int count;
    if (myinstream.good()) {
      (myinstream)>>count;
      std::string buffer;
	(myinstream)>>m_id;
	(myinstream)>>iy>>idum2;
	for (int i=0;i<NTAB;++i) (myinstream)>>iv[i];
      myinstream.close();
    } 
    else msg_Error()<<"ERROR in Random::ReadInStatus : "<<filename<<" not found!!"<<endl;
  }
}

size_t ATOOLS::Random::ReadInStatus
(std::istream &myinstream,const size_t &idx)
{
  if (activeGenerator==4) { ReadInStatus4(myinstream,idx); } 
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_size()>1) return std::string::npos;
#endif
  size_t count;
  while (!myinstream.eof()) {
    myinstream>>count;
    std::string buffer;
    myinstream>>m_id;
    myinstream>>iy>>idum2;
    for (int i=0;i<NTAB;++i) myinstream>>iv[i];
    if (count==idx) return idx+1;
  }
  return std::string::npos;
}

double ATOOLS::Random::GetNZ() 
{
  double ran1;
  do ran1=Get(); while (ran1==0.); 
  return ran1;
}


void ATOOLS::Random::SetSeed(long int nid) 
{
  m_id = nid<0 ? nid : -nid;
  activeGenerator = 2;
}


void ATOOLS::Random::SaveStatus()
{
  if (p_external!=NULL) return;
  if (activeGenerator==4) { return SaveStatus4(); };
  m_sid=m_id; 
  siy=iy;
  sidum2=idum2;
  for (int i=0;i<NTAB;++i) siv[i]=iv[i];
}


void ATOOLS::Random::RestoreStatus()
{
  if (p_external!=NULL) return;
  if (activeGenerator==4) { return RestoreStatus4(); };
  m_id=m_sid; 
  iy=siy;
  idum2=sidum2;
  for (int i=0;i<NTAB;++i) iv[i]=siv[i];
}

void ATOOLS::Random:: ResetToLastIncrementedSeed()
{
  m_nsinceinit=0;
  ReadInStatus(m_lastincrementedseed,0);
}


void ATOOLS::Random::PrepareTerminate()
{
  if (p_external!=NULL) return;
  std::string path(rpa->gen.Variable("SHERPA_STATUS_PATH"));
  if (path=="") return;
  RestoreStatus();
  WriteOutStatus((path+"/random.dat").c_str());
}

// ----------------- Methods for new Random Number Generator -------------

/*
   This is the random number generator proposed by George Marsaglia in
   Florida State University Report: FSU-SCRI-87-50
*/
/*
   This is the initialization routine for the random number generator.
   NOTE: The seed variables can have values between:    0 <= IJ <= 31328
                                                        0 <= KL <= 30081
   The random number sequences created by these two seeds are of sufficient
   length to complete an entire calculation with. For example, if sveral
   different groups are working on different parts of the same calculation,
   each group could be assigned its own IJ seed. This would leave each group
   with 30000 choices for the second seed. That is to say, this random
   number generator can create 900 million different subsequences -- with
   each subsequence having a length of approximately 10^30.
*/

bool Random::InitExternal(const std::string &path,const std::string &file)
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.SetInputPath(path);
  read.SetInputFile(file);
  std::string name;
  if (!read.ReadFromFile(name,"EXTERNAL_RNG")) return false;
  p_external = RNG_Getter::GetObject(name,RNG_Key());
  if (p_external==NULL) {
    msg_Out()<<METHOD<<"(): {\n\n  // available RNGs\n\n";
    RNG_Getter::PrintGetterInfo(msg->Out(),20);
    msg_Out()<<"\n\n}"<<std::endl;
    THROW(fatal_error,"RNG '"+name+"' not found.");
  }
  msg_Info()<<METHOD<<"(): Initialized external RNG '"<<name<<"'."<<std::endl;
  return true;
}

ATOOLS::Random::Random(unsigned int i1,unsigned int i2,unsigned int i3,
		       unsigned int i4,unsigned int i5,unsigned int i6) : 
  m_lastincrementedseed(ios_base::binary|ios_base::in|ios_base::out),
  m_nsinceinit(0), m_increment(0), p_external(NULL)
{  
  ATOOLS::exh->AddTerminatorObject(this);
  p_ran4[0] = new Marsaglia();
  SetSeed(i1,i2,i3,i4);
  p_ran4[1] = new Marsaglia(*p_ran4[0]);
}


void ATOOLS::Random::SetSeed(unsigned int i1,unsigned int i2,
			     unsigned int i3,unsigned int i4)
{
  msg_Info()<<METHOD<<"(): Seeds set to "<<i1<<" "<<i2
	    <<" "<<i3<<" "<<i4<<std::endl;
  p_ran4[0]->Init(i1,i2,i3,i4);
  *p_ran4[1]=*p_ran4[0];
  activeGenerator = 4;
}


double ATOOLS::Random::Ran4()
{
  return p_ran4[0]->Get();
}


int ATOOLS::Random::WriteOutStatus4(const char * filename)
{
  // remove old random file
  if (FileExists(filename)) remove(filename);
  // open file and append status of Random Generator at the end if possible
  std::ofstream file((std::string(filename)+".msg").c_str());
  p_ran4[0]->WriteStatus(file);
  return 1;
}

int ATOOLS::Random::WriteOutSavedStatus4(const char * filename)
{
  // remove old random file
  if (FileExists(filename)) remove(filename);
  // open file and append status of Random Generator at the end if possible
  std::ofstream file((std::string(filename)+".msg").c_str());
  p_ran4[1]->WriteStatus(file);
  return 1;
}

int ATOOLS::Random::WriteOutStatus4(std::ostream &os,const size_t &idx)
{
  THROW(fatal_error,"Invalid call");
  return -1;
}

void ATOOLS::Random::ReadInStatus4(const char * filename)
{
  msg_Info()<<"Random::ReadInStatus from "<<filename<<".msg"<<endl;

  std::ifstream file((std::string(filename)+".msg").c_str());
  if (file.good()) p_ran4[0]->ReadStatus(file);  
  else 
    msg_Error()<<"ERROR in Random::ReadInStatus4 : "<<filename<<" not found!!"<<endl;
  *p_ran4[1]=*p_ran4[0];
}

size_t ATOOLS::Random::ReadInStatus4(std::istream &is,const size_t &idx)
{
  THROW(fatal_error,"Invalid call");
  return 0;
}

void ATOOLS::Random::SaveStatus4() { *p_ran4[1]=*p_ran4[0]; }
void ATOOLS::Random::RestoreStatus4() { *p_ran4[0]=*p_ran4[1]; }

double ATOOLS::Random::Get()   
{
  if (p_external) return p_external->Get();
  // Sherpa internal
  double rng(0.);
  ++m_nsinceinit;
  if (activeGenerator==4) rng=Ran4();
  else                    rng=Ran2(&m_id);
  if (activeGenerator!=4 && m_increment && m_nsinceinit==m_increment) {
    EraseLastIncrementedSeed();
    WriteOutStatus(m_lastincrementedseed,0);
  }
  return rng;
}

ptrdiff_t ATOOLS::Random::operator() (ptrdiff_t max) {
  return Min(static_cast<ptrdiff_t>(Get() * max),max-1);
}

void ATOOLS::Random::Gaussian(double & x,double & y)   
{
  double phi(2.*M_PI*Get()), random(Get());
  while (random==0.) random = Get();
  double r(sqrt(-2.*log(random)));

  x = r*std::cos(phi);
  y = r*std::sin(phi);
}

External_RNG::~External_RNG()
{
}
