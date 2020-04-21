#include "PHASIC++/Channels/Exponential_Channels.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Exponential_RelicDensity::
Exponential_RelicDensity(const double exponent,const double mass1,const double mass2,
        const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_exponent(exponent)
{
  m_mass[0]=mass1;
  m_mass[1]=mass2;
  m_name="Exponential_"+ATOOLS::ToString(exponent)+"_RelicDensity";
  m_spkey.SetInfo(std::string("Exponential_")+ATOOLS::ToString(exponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_zchannel = m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum   = 1;
  p_vegas    = new Vegas(m_rannum,100,m_name,0);
  p_rans     = new double[m_rannum];
}

void Exponential_RelicDensity::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  m_spkey[3] = CE.ExponentialMomenta(m_exponent,m_spkey[0],m_spkey[1],m_mass,p_rans[0]);
}

void Exponential_RelicDensity::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ExponentialWeight(m_exponent,m_spkey[0],m_spkey[1],m_mass,
          m_spkey[3],m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  p_rans[0] = m_sgridkey[0];
  double pw = p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight();
}

///////////////////////////////////////////////////////////////////////////////

Exponential_DM_Annihilation::
Exponential_DM_Annihilation(const double exponent,const double mass1,const double mass2,
        const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_exponent(exponent)
{
  m_mass[0]=mass1;
  m_mass[1]=mass2;
  m_name="Exponential_"+ATOOLS::ToString(exponent)+"_DM_Annihilation";
  m_spkey.SetInfo(std::string("Exponential_")+ATOOLS::ToString(exponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_xkey.Assign(cinfo+std::string("::xDM"),3,0,info);
  m_cosxikey.Assign(cinfo+std::string("::cosXi"),3,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_xgridkey.Assign(m_xkey.Info(),1,0,info);
  m_cosgridkey.Assign(m_cosxikey.Info(),1,0,info);
  m_zchannel = m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum   = 2;
  p_vegas    = new Vegas(m_rannum,100,m_name,0);
  p_rans     = new double[m_rannum];
}

void Exponential_DM_Annihilation::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  m_spkey[3] = m_spkey[2] = CE.ExponentialMomenta(m_exponent,m_spkey[0],m_spkey[1],m_mass,p_rans[0]);
  m_cosxikey[2] = CE.GenerateDMAngleUniform(p_rans[1],3);
  m_xkey[2] = CE.GenerateDMRapidityUniform(m_mass,m_spkey.Doubles(),m_xkey.Doubles(),
					   m_cosxikey[2], p_rans[0], 3);
}

void Exponential_DM_Annihilation::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ExponentialWeight(m_exponent,m_spkey[0],m_spkey[1],m_mass,
          m_spkey[3],m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  p_rans[0] = m_sgridkey[0];
  double pw = p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()/m_spkey[2];
  // msg_Out()<<"s="<<m_spkey[3]<<", weight="<<m_weight<<"\n"; //debugging
}
