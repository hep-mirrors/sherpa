#include "ATOOLS/Math/Marsaglia.H"

#include <cstdlib>

using namespace ATOOLS;

#define UC (unsigned char)
#define znew(z) (z=36969*(z&65535)+(z>>16))
#define wnew(w) (w=18000*(w&65535)+(w>>16))
#define MWC(z,w) ((znew(z)<<16)+wnew(w))
#define SHR3(jsr) (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define CONG(jcong) (jcong=69069*jcong+1234567)
#define KISS(z,w,jcong,jsr) ((MWC(z,w)^CONG(jcong))+SHR3(jsr))
#define SWB(c,bro,x,y,t) (c++,bro=(x<y),t[c]=(x=t[UC(c+34)])-(y=t[UC(c+19)]+bro))

Marsaglia::Marsaglia(): m_x(0), m_y(0), m_c(0)
{
  if (sizeof(UL)!=4) {
    std::cout<<"Unsigned int type has size "
	     <<sizeof(UL)<<", should be 4."<<std::endl;
    exit(1);
  }
  Init(12345,65435,34221,12345);
  for(int i=1;i<1000000;++i) SWB(m_c,m_bro,m_x,m_y,m_t);
  if (SWB(m_c,m_bro,m_x,m_y,m_t)!=1429146441U) {
    std::cout<<"RNG test 1 failed."<<std::endl;
    exit(1);
  }
  for(int i=1;i<1000000;++i) KISS(m_z,m_w,m_jcong,m_jsr);
  if (KISS(m_z,m_w,m_jcong,m_jsr)!=1372460312U) {
    std::cout<<"RNG test 2 failed."<<std::endl;
    exit(1);
  }
}

void Marsaglia::Init(UL i1,UL i2,UL i3,UL i4)
{
  m_z=i1, m_w=i2, m_jsr=i3, m_jcong=i4;
  for(int i=0;i<256;++i) m_t[i]=KISS(m_z,m_w,m_jcong,m_jsr);
}

double Marsaglia::Get()
{
  // original: (KISS+SWB)/4294967295
  // -> draws in interval [0,1] including exact 0. and 1.
  // change to (KISS+SWB+1)/4294967297
  // -> draws in interval [1/4294967297,4294967296/4294967297]
  return (KISS(m_z,m_w,m_jcong,m_jsr)
          +SWB(m_c,m_bro,m_x,m_y,m_t)+1.)*2.328306436e-10;
}

void Marsaglia::WriteStatus(std::ostream &str)
{
  str.write((char*)this,sizeof(*this));
}

bool Marsaglia::ReadStatus(std::istream &str)
{
  str.read((char*)this,sizeof(*this));
  return true;
}
