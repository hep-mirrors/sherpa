#include "PHASIC++/Channels/Threshold.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <stdio.h>

using namespace PHASIC;

Threshold_Uniform::Threshold_Uniform(const double mass,const std::string cinfo,
				     ATOOLS::Integration_Info *info):
  m_mass(mass)
{
  name=std::string("Threshold_Uniform_")+ATOOLS::ToString((int)(100.*mass));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Threshold_Uniform::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rans,const int mode) 
{
  m_spkey[3]=CE.ThresholdMomenta(m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Uniform::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_mass,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

Threshold_Forward::Threshold_Forward(const double mass,const double yexponent,
				     const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_mass(mass), 
  m_yexponent(yexponent)
{
  name=std::string("Threshold_Forward_")+ATOOLS::ToString((int)(100.*mass));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Threshold_Forward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rans,const int mode) 
{
  m_spkey[3]=CE.ThresholdMomenta(m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Forward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_mass,m_spkey[0],m_spkey[1],m_spkey[3]);
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

Threshold_Backward::Threshold_Backward(const double mass,const double yexponent,
				       const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_mass(mass), 
  m_yexponent(yexponent)
{
  name=std::string("Threshold_Backward_")+ATOOLS::ToString((int)(100.*mass));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Threshold_Backward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				       const double *rans,int mode)
{
  m_spkey[3]=CE.ThresholdMomenta(m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Backward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_mass,m_spkey[0],m_spkey[1],m_spkey[3]);
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

Threshold_Central::Threshold_Central(const double mass,const std::string cinfo,
				     ATOOLS::Integration_Info *info):
  m_mass(mass)
{
  name=std::string("Threshold_Central_")+ATOOLS::ToString((int)(100.*mass));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Threshold_Central::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rans,int mode)
{
  m_spkey[3]=CE.ThresholdMomenta(m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Central::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_mass,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

