#include "LL_KPerp.H"

#include "Channel_Elements.H"

#include <stdio.h>

using namespace PHASIC;

LL_KPerp::LL_KPerp(double _m_beta,const std::string cinfo,
		   Integration_Info *info):
  m_beta(_m_beta)
{
  char beta[3];
  sprintf(beta,"%i",int(m_beta*100.));
  name=std::string("LL_KPerp_")+std::string(beta);
  m_kp1key.SetInfo(name+cinfo);
  m_kp2key.SetInfo(name+cinfo);
  m_kp1key.Assign("k_perp_1",3,0,info);
  m_kp2key.Assign("k_perp_2",3,0,info);
}

void LL_KPerp::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode)
{
  CalculateLimits(spkey,ykey);
  m_kp1key[2]=CE.MasslessPropMomenta(m_beta,m_kp1key[0],m_kp1key[1],rans[0]);
  m_kp2key[2]=CE.MasslessPropMomenta(m_beta,m_kp2key[0],m_kp2key[1],rans[1]);
}

void LL_KPerp::GenerateWeight(const int mode)
{
  if (m_kp1key[2]>=m_kp1key[0] && m_kp1key[2]<=m_kp1key[1]) {
    if (m_kp1key.Weight()==UNDEFINED_WEIGHT) {
      m_kp1key<<1./CE.MasslessPropWeight(m_beta,m_kp1key[0],m_kp1key[1],m_kp1key[2]);
    }
  }
  if (m_kp2key[2]>=m_kp2key[0] && m_kp2key[2]<=m_kp2key[1]) {
    if (m_kp2key.Weight()==UNDEFINED_WEIGHT) {
      m_kp2key<<1./CE.MasslessPropWeight(m_beta,m_kp2key[0],m_kp2key[1],m_kp2key[2]);
    }
  }
  weight=m_kp1key.Weight()*m_kp2key.Weight();
}

void LL_KPerp::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  m_kp1key[1]=spkey[3];
  m_kp2key[1]=spkey[3];
}

