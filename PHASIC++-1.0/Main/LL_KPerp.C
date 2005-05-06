#include "LL_KPerp.H"

#include "Channel_Elements.H"
#include "MyStrStream.H"
#include "Random.H"

#include <stdio.h>

using namespace PHASIC;

LL_KPerp::LL_KPerp(double beta,const std::string cinfo,
		   ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_calculated(false)
{
  name=std::string("LL_KPerp_")+ATOOLS::ToString(beta);
  m_kpkey[0].SetInfo(name+cinfo);
  m_kpkey[1].SetInfo(name+cinfo);
  m_kpkey[0].Assign("k_perp_1",4,1,info);
  m_kpkey[1].Assign("k_perp_2",4,1,info);
  rannum=2;
  p_vegas = new Vegas(2,250,name,0);
  rans = new double[2];
}

LL_KPerp::~LL_KPerp()
{
  delete [] rans;
  delete p_vegas;
}

double LL_KPerp::f(ATOOLS::Info_Key &key,double t)
{
  if (t<key[1]) return pow(key[1],-m_beta);
  return pow(t,-m_beta);
}

double LL_KPerp::I(ATOOLS::Info_Key &key,double t)
{
  if (t<key[1]) return (t-key[0])*pow(key[1],-m_beta);
  double I1=(key[1]-key[0])*pow(key[1],-m_beta);
  return I1+(pow(t,1.-m_beta)-pow(key[1],1.-m_beta))/(1.-m_beta);
}

double LL_KPerp::t(ATOOLS::Info_Key &key,double I)
{
  double I1=(key[1]-key[0])*pow(key[1],-m_beta);
  if (I<I1) return I/pow(key[1],-m_beta)+key[0];
  I-=I1;
  return pow(I*(1.-m_beta)+pow(key[1],1.-m_beta),1.0/(1.-m_beta));
}

std::string LL_KPerp::ChID()
{
  return name;
}

void LL_KPerp::WriteOut(std::string pId) 
{ 
  p_vegas->WriteOut(pId); 
}

void LL_KPerp::ReadIn(std::string pId)   
{ 
  p_vegas->ReadIn(pId); 
}

void LL_KPerp::Optimize()  
{
  p_vegas->Optimize();
} 

void LL_KPerp::EndOptimize()  
{
  p_vegas->EndOptimize();
} 

void LL_KPerp::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
			     const double *rans,const int mode)
{
  CalculateLimits(spkey,ykey);
  for (short unsigned int i=0;i<2;++i) {
    m_kpkey[i][3]=t(m_kpkey[i],I(m_kpkey[i],m_kpkey[i][0])*(1.0-rans[i])+
		    rans[i]*I(m_kpkey[i],m_kpkey[i][2]));
    double phi=2.0*M_PI*ATOOLS::ran.Get();
    double kp=sqrt(m_kpkey[i][3]); 
    m_kpkey[i](0)=ATOOLS::Vec4D(0.0,cos(phi)*kp,sin(phi)*kp,0.0);
  }
}

void LL_KPerp::GenerateWeight(const int mode)
{
  for (short unsigned int i=0;i<2;++i) 
    if (m_kpkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) 
      m_kpkey[i]<<(I(m_kpkey[i],m_kpkey[i][2])-I(m_kpkey[i],m_kpkey[i][0]))/
	f(m_kpkey[i],m_kpkey[i][3]);
  weight=m_kpkey[0].Weight()*m_kpkey[1].Weight();
}

void LL_KPerp::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  for (short unsigned int i=0;i<2;++i) 
    rans[i]=(I(m_kpkey[i],m_kpkey[i][3])-I(m_kpkey[i],m_kpkey[i][0]))/
      (I(m_kpkey[i],m_kpkey[i][2])-I(m_kpkey[i],m_kpkey[i][0]));
  p_vegas->AddPoint(value,rans);
}

void LL_KPerp::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
  for (short unsigned int i=0;i<2;++i) m_kpkey[1][2]=spkey[2];
}

