#include "PHASIC++/Channels/LBS_Compton_Peak_Channels.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

LBS_Compton_Peak_Uniform::
LBS_Compton_Peak_Uniform(const double exponent,const double pole,
			 const std::string cinfo,ATOOLS::Integration_Info *info,
			 const size_t mode):
  ISR_Channel_Base(info),
  m_exponent(exponent),
  m_pole(pole), m_mode(mode)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  m_name=std::string("LBS_Compton_Peak_Uniform");
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel = m_spkey.Name().find("z-channel")!=std::string::npos;
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_rannum = 2;
  p_vegas  = new Vegas(2,100,m_name);
  p_rans   = new double[2];
}

void LBS_Compton_Peak_Uniform::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],p_rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYUniform(sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),p_rans[1],m_mode);
}

void LBS_Compton_Peak_Uniform::GenerateWeight(const int & mode)
{
  m_weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,
				  m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYUniform(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),
				m_ygridkey[0],m_mode);
    }
  }
  p_rans[0] = m_sgridkey[0];
  p_rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



LBS_Compton_Peak_Forward::
LBS_Compton_Peak_Forward(const double exponent,const double pole,
			 const double yexponent,
			 const std::string cinfo,ATOOLS::Integration_Info *info,
			 const size_t mode):
  ISR_Channel_Base(info),
  m_exponent(exponent),
  m_pole(pole),
  m_yexponent(yexponent), m_mode(mode)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  m_name=std::string("LBS_Compton_Peak_Forward");
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_rannum=2;
  p_vegas = new Vegas(2,100,m_name);
  p_rans  = new double[2];
}

void LBS_Compton_Peak_Forward::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],p_rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYForward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),p_rans[1],m_mode);
}

void LBS_Compton_Peak_Forward::GenerateWeight(const int & mode)
{
  m_weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,
				  m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYForward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],m_mode);
    }
  }
  p_rans[0] = m_sgridkey[0];
  p_rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


LBS_Compton_Peak_Backward::
LBS_Compton_Peak_Backward(const double exponent,const double pole,
			  const double yexponent,
			  const std::string cinfo,ATOOLS::Integration_Info *info,
			  const size_t mode):
  ISR_Channel_Base(info),
  m_exponent(exponent),
  m_pole(pole),
  m_yexponent(yexponent), m_mode(mode)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  m_name=std::string("LBS_Compton_Peak_Backward");
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_rannum=2;
  p_vegas = new Vegas(2,100,m_name);
  p_rans  = new double[2];
}

void LBS_Compton_Peak_Backward::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],p_rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYBackward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),p_rans[1],m_mode);
}

void LBS_Compton_Peak_Backward::GenerateWeight(const int & mode)
{
  m_weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,
				  m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYBackward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],m_mode);
    }
  }
  p_rans[0] = m_sgridkey[0];
  p_rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


LBS_Compton_Peak_Central::
LBS_Compton_Peak_Central(const double exponent,const double pole,
			 const std::string cinfo,ATOOLS::Integration_Info *info,
			 const size_t mode):
  ISR_Channel_Base(info),
  m_exponent(exponent),
  m_pole(pole), m_mode(mode)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  m_name=std::string("LBS_Compton_Peak_Central");
  m_ykey.SetInfo("Central");
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_rannum=2;
  p_vegas = new Vegas(m_rannum,100,m_name);
  p_rans  = new double[2];
}

void LBS_Compton_Peak_Central::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],p_rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYCentral(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),
				p_rans[1],m_mode);
}

void LBS_Compton_Peak_Central::GenerateWeight(const int & mode)
{
  m_weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,
				  m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYCentral(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),
				m_ygridkey[0],m_mode);
    }
  }
  p_rans[0] = m_sgridkey[0];
  p_rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

