#include "LL_KPerp.H"

#include "Channel_Elements.H"
#include "Scaling.H"

#include <stdio.h>

using namespace PHASIC;

LL_KPerp::LL_KPerp(double beta,const std::string cinfo,
		   ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_calculated(false)
{
  name=std::string("LL_KPerp_")+ATOOLS::ToString(beta);
  m_kp1key.SetInfo(name+cinfo);
  m_kp2key.SetInfo(name+cinfo);
  m_kp1key.Assign("k_perp_1",4,0,info);
  m_kp2key.Assign("k_perp_2",4,0,info);
}

void LL_KPerp::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode)
{
  CalculateLimits(spkey,ykey);
  m_kp1key[3]=rans[0]*m_integral[0];
  if (rans[0]>m_rancut[0]) {
    m_kp1key[3]=m_kp1key[3]-m_firstintegral[0]; 
    m_kp1key[3]=pow(m_kp1key[3]*(1.-m_beta)+pow(m_kp1key[1],1.-m_beta),1./(1.-m_beta));
  }
  else {
    m_kp1key[3]=m_kp1key[3]*pow(m_kp1key[1],m_beta)+m_kp1key[0];
  }
  m_kp2key[3]=rans[1]*m_integral[1];
  if (rans[1]>m_rancut[1]) {
    m_kp2key[3]=m_kp2key[3]-m_firstintegral[1];
    m_kp2key[3]=pow(m_kp2key[3]*(1.-m_beta)+pow(m_kp2key[1],1.-m_beta),1./(1.-m_beta));
  }
  else {
    m_kp2key[3]=m_kp2key[3]*pow(m_kp2key[1],m_beta)+m_kp2key[0];
  }
}

void LL_KPerp::GenerateWeight(const int mode)
{
  CalculateLimits();
  if (m_kp1key[3]>=m_kp1key[0] && m_kp1key[3]<=m_kp1key[2]) {
    if (m_kp1key.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double dist=1.;
      if (m_kp1key[3]>m_kp1key[1]) dist=pow(m_kp1key[3],m_beta);
      else dist=pow(m_kp1key[1],m_beta);
      m_kp1key<<m_integral[0]*dist;
    }
  }
  if (m_kp2key[3]>=m_kp2key[0] && m_kp2key[3]<=m_kp2key[2]) {
    if (m_kp2key.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double dist=1.;
      if (m_kp2key[3]>m_kp2key[1]) dist=pow(m_kp2key[3],m_beta);
      else dist=pow(m_kp2key[1],m_beta);
      m_kp2key<<m_integral[1]*dist;
    }
  }
  weight=m_kp1key.Weight()*m_kp2key.Weight();
  // m_calculated=false;
}

void LL_KPerp::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
  m_kp1key[2]=spkey[2];
  m_kp2key[2]=spkey[2];
  CalculateLimits();
}

void LL_KPerp::CalculateLimits() 
{
  if (m_calculated) return;
  m_firstintegral[0]=(m_kp1key[1]-m_kp1key[0])/pow(m_kp1key[1],m_beta);
  m_integral[0]=m_firstintegral[0]+(pow(m_kp1key[2],1.-m_beta)-pow(m_kp1key[1],1.-m_beta))/(1.-m_beta);
  m_firstintegral[1]=(m_kp2key[1]-m_kp2key[0])/pow(m_kp2key[1],m_beta);
  m_integral[1]=m_firstintegral[1]+(pow(m_kp2key[2],1.-m_beta)-pow(m_kp2key[1],1.-m_beta))/(1.-m_beta);
  m_rancut[0]=m_firstintegral[0]/m_integral[0];
  m_rancut[1]=m_firstintegral[1]/m_integral[1];
  m_calculated=true;
}

