#include "Random.H"
#include "Message.H"
#include "MathTools.H"
#include <iostream>

using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

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

AMATOOLS::Random AMATOOLS::Ran(1234);

double Random::ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
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
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double Random::ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
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
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

void Random::Init_Ran3(long *idum)
{
  long mj,mk;
  int i,ii,k;
  
  mj=MSEED-(*idum < 0 ? -*idum : *idum);
  mj %= MBIG;
  ma[55]=mj;
  mk=1;
  for (i=1;i<=54;i++) {
    ii=(21*i) % 55;
    ma[ii]=mk;
    mk=mj-mk;
    if (mk < MZ) mk += MBIG;
    mj=ma[ii];
  }
  for (k=1;k<=4;k++)
    for (i=1;i<=55;i++) {
      ma[i] -= ma[1+(i+30) % 55];
      if (ma[i] < MZ) ma[i] += MBIG;
    }
  inext=0;
  inextp=31;
}

double Random::Ran3()
{
  //  cout<<" inext= "<<inext<<endl;
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  long mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */


int Random::WriteOutStatus(char * filename){
  // write out every Statusregister of Random Number generator

  //  sprintf(outname,"%s%i.dat",filename,written); 
  if ((myoutstream!=0) && (strcmp(filename,outname)!=0)) {
    myoutstream->close();
    myoutstream = 0;
  }
  if (myoutstream == 0){
    msg.Events()<<" Saving Random Number Generator Status to "<<filename<<endl;
    long int count=0;
    std::ifstream myinstream(filename);
    if (!myinstream.bad()) {
      char * buffer[600];
      while (!myinstream.eof()) {
	myinstream>>count;
	myinstream.getline((char*)buffer,600);
      }
      ++count;
      myinstream.close();
    }
#ifdef _IOS_BAD
    myoutstream=new std::fstream(filename,ios::app);
#else
    myoutstream=new std::fstream(filename,std::ios_base::app);
#endif
    strcpy(outname,filename);
    written=count;
  } 
  (*myoutstream)<<written<<"\t"<<id<<"\t"<<inext<<"\t"<<inextp<<"\t";
  for (int i=0;i<56;++i) (*myoutstream)<<ma[i]<<"\t";
  (*myoutstream)<<endl;
  return written++;
}

void Random::ReadInStatus(char * filename, long int index){
  // read in every Statusregister of Random Number generator
  msg.Events()<<"reading file "<<filename<<" index "<<index<<endl;
  std::ifstream myinstream(filename);
  long int count;
  if (!myinstream.bad()) {
    (myinstream)>>count;
    char buffer[600];
    while ((count!=index)&&(!myinstream.eof())) {
      myinstream.getline(buffer,600);
      (myinstream)>>count;    
    };
    if (count==index) {
      msg.Events()<<" index="<<index<<" found"<<endl;
      (myinstream)>>id; (myinstream)>>inext; (myinstream)>>inextp;
      for (int i=0;i<56;++i) (myinstream)>>ma[i];    
    } 
    else msg.Error()<<" index="<<index<<" not found in "<<filename<<endl;
    myinstream.close();
  } 
  else msg.Error()<<filename<<" not found!!"<<endl;
}
