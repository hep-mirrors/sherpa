#include "LDL_KPerp.H"

#include <stdio.h>

using namespace PHASIC;

LDL_KPerp::LDL_KPerp(const double lambda2,const std::string cinfo,
		     ATOOLS::Integration_Info *info):
  m_lambda2(lambda2)
{
  char lambda[3];
  sprintf(lambda,"%i",int(m_lambda2*10000.));
  name=std::string("LDL_KPerp_")+std::string(lambda);
  m_kp1key.SetInfo(name+cinfo);
  m_kp2key.SetInfo(name+cinfo);
  m_kp1key.Assign("k_perp_1",4,0,info);
  m_kp2key.Assign("k_perp_2",4,0,info);
}

void LDL_KPerp::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode)
{
  CalculateLimits(spkey,ykey);
  if (rans[0]>m_rancut[0]) {
    m_kp1key[3]=m_lambda2*pow(m_kp1key[1]/m_lambda2,
			      pow(log(m_kp1key[2]/m_lambda2)/log(m_kp1key[1]/m_lambda2),
				  (rans[0]-m_rancut[0])/(1.-m_rancut[0])));
  }
  else {
    m_kp1key[3]=m_kp1key[0]+(m_kp1key[1]-m_kp1key[0])*rans[0]/m_rancut[0];
  }
  if (rans[1]>m_rancut[1]) {
    m_kp2key[3]=m_lambda2*pow(m_kp2key[1]/m_lambda2,
			      pow(log(m_kp2key[2]/m_lambda2)/log(m_kp2key[1]/m_lambda2),
				  (rans[1]-m_rancut[1])/(1.-m_rancut[1])));
  }
  else {
    m_kp2key[3]=m_kp2key[0]+(m_kp2key[1]-m_kp2key[0])*rans[1]/m_rancut[1];
  }
}

void LDL_KPerp::GenerateWeight(const int mode)
{
  if (m_kp1key[3]>=m_kp1key[0] && m_kp1key[3]<=m_kp1key[2]) {
    if (m_kp1key.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      if (m_kp1key[3]>m_kp1key[1]) {
	m_kp1key<<log(log(m_kp1key[2]/m_lambda2)/log(m_kp1key[1]/m_lambda2))*
	  m_kp1key[3]*log(m_kp1key[3]/m_lambda2);
      }
      else {
	m_kp1key<<m_kp1key[1]-m_kp1key[0];
      }
    }
  }
  if (m_kp2key[3]>=m_kp2key[0] && m_kp2key[3]<=m_kp2key[2]) {
    if (m_kp2key.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      if (m_kp2key[3]>m_kp2key[1]) {
	m_kp2key<<log(log(m_kp2key[2]/m_lambda2)/log(m_kp2key[1]/m_lambda2))*
	  m_kp2key[3]*log(m_kp2key[3]/m_lambda2);
      }
      else {
	m_kp2key<<m_kp2key[1]-m_kp2key[0];
      }
    }
  }
  weight=m_kp1key.Weight()*m_kp2key.Weight();
}

void LDL_KPerp::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
  m_kp1key[2]=spkey[2];
  m_kp2key[2]=spkey[2];
  m_rancut[0]=(m_kp1key[1]-m_kp1key[0])/(m_kp1key[2]-m_kp1key[0]);
  m_rancut[1]=(m_kp2key[1]-m_kp2key[0])/(m_kp2key[2]-m_kp2key[0]);
}

