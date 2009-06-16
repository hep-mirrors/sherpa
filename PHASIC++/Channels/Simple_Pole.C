#include "PHASIC++/Channels/Simple_Pole.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <stdio.h>

using namespace PHASIC;

Simple_Pole_Uniform::Simple_Pole_Uniform(const double exponent,const std::string cinfo,
					 ATOOLS::Integration_Info *info):
  m_exponent(exponent)
{
  name=std::string("Simple_Pole_Uniform_")+ATOOLS::ToString((int)(100.*exponent));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Simple_Pole_Uniform::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rans,const int mode) 
{
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Uniform::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

Simple_Pole_Forward::Simple_Pole_Forward(const double sexponent,const double yexponent,
					 const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_sexponent(sexponent), 
  m_yexponent(yexponent)
{
  name=std::string("Simple_Pole_Forward_")+ATOOLS::ToString((int)(100.*sexponent));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(sexponent));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Simple_Pole_Forward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rans,const int mode) 
{
  m_spkey[3]=CE.MasslessPropMomenta(m_sexponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Forward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_sexponent,m_spkey[0],m_spkey[1],m_spkey[3]);
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

Simple_Pole_Backward::Simple_Pole_Backward(const double sexponent,const double yexponent,
					   const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_sexponent(sexponent), 
  m_yexponent(yexponent)
{
  name=std::string("Simple_Pole_Backward_")+ATOOLS::ToString((int)(100.*sexponent));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(sexponent));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Simple_Pole_Backward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					 const double *rans,int mode)
{
  m_spkey[3]=CE.MasslessPropMomenta(m_sexponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Backward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_sexponent,m_spkey[0],m_spkey[1],m_spkey[3]);
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

Simple_Pole_Central::Simple_Pole_Central(const double exponent,const std::string cinfo,
					 ATOOLS::Integration_Info *info):
  m_exponent(exponent)
{
  name=std::string("Simple_Pole_Central_")+ATOOLS::ToString((int)(100.*exponent));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Simple_Pole_Central::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rans,int mode)
{
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Central::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}


