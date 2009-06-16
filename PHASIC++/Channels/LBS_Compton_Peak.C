#include "PHASIC++/Channels/LBS_Compton_Peak.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <stdio.h>

using namespace PHASIC;

LBS_Compton_Peak_Uniform::LBS_Compton_Peak_Uniform(const double exponent,const double pole,
						   const std::string cinfo,ATOOLS::Integration_Info *info):
  m_exponent(exponent),
  m_pole(pole)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  name=std::string("LBS_Compton_Peak_Uniform");
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void LBS_Compton_Peak_Uniform::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					     const double *rans,const int mode) 
{
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Uniform::GenerateWeight(const int mode) 
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      double help=m_spkey[3];
      if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
	if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
	else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
      }
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

LBS_Compton_Peak_Forward::LBS_Compton_Peak_Forward(const double exponent,const double pole,
						   const double yexponent,
						   const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_exponent(exponent),
  m_pole(pole),
  m_yexponent(yexponent)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  name=std::string("LBS_Compton_Peak_Forward");
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void LBS_Compton_Peak_Forward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					     const double *rans,const int mode) 
{
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Forward::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      double help=m_spkey[3];
      if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
	if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
	else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
      }
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

LBS_Compton_Peak_Backward::LBS_Compton_Peak_Backward(const double exponent,const double pole,
						     const double yexponent,
						     const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_exponent(exponent),
  m_pole(pole),
  m_yexponent(yexponent)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  name=std::string("LBS_Compton_Peak_Backward");
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void LBS_Compton_Peak_Backward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					      const double *rans,int mode)
{
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Backward::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      double help=m_spkey[3];
      if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
	if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
	else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
      }
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

LBS_Compton_Peak_Central::LBS_Compton_Peak_Central(const double exponent,const double pole,
						   const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_exponent(exponent),
  m_pole(pole)
{
  std::string help=ATOOLS::ToString(exponent)+
    std::string("_")+ATOOLS::ToString(pole);
  m_spkey.SetInfo(std::string("LBS_Compton_Peak_")+help);
  name=std::string("LBS_Compton_Peak_Central");
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void LBS_Compton_Peak_Central::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					     const double *rans,int mode)
{
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Central::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      double help=m_spkey[3];
      if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
	if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
	else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
      }
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

