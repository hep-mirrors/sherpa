#include "Resonance.H"

#include "Channel_Elements.H"
#include "Message.H"

#include <stdio.h>

using namespace PHASIC;

Resonance_Uniform::Resonance_Uniform(const double mass,const double width,
				     const std::string cinfo,Integration_Info *info):
  m_mass(mass),
  m_width(width)
{
  char help[3];
  sprintf(help,"%i",int(mass*1000.));
  name=std::string(help);
  sprintf(help,"%i",int(width*1000.));
  name+=std::string(help);
  m_spkey.SetInfo(std::string("Resonance_")+std::string(help));
  name=std::string("Resonance_Uniform");
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Resonance_Uniform::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]+=CE.DiceYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Resonance_Uniform::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

bool Resonance_Uniform::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  m_spkey[2]=spkey[2];
  if (!m_zchannel) {
    m_spkey[0]=spkey[0];
    m_spkey[1]=spkey[3];
  }
  else {
    m_spkey[0]=spkey[3];
    m_spkey[1]=spkey[1];
    m_ykey[0]=ykey[0];
    m_ykey[1]=ykey[1];
    double logtau=.5*log(spkey[3]/spkey[2]);
    m_xkey[0]=logtau+ykey[2];
    m_xkey[2]=logtau-ykey[2];
    m_xkey[1]=0.;
    m_xkey[3]=0.;
  }
  return true;
}

Resonance_Forward::Resonance_Forward(const double mass,const double width,const double yexponent,
				   const std::string cinfo,Integration_Info *info): 
  m_mass(mass),
  m_width(width),
  m_yexponent(yexponent)
{
  char help[3];
  sprintf(help,"%i",int(mass*1000.));
  name=std::string(help);
  sprintf(help,"%i",int(width*1000.));
  name+=std::string(help);
  m_spkey.SetInfo(std::string("Resonance_")+std::string(help));
  name=std::string("Resonance_Forward");
  sprintf(help,"%i",int(100.*yexponent));
  m_ykey.SetInfo(std::string("Forward_")+std::string(help));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Resonance_Forward::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]+=CE.DiceYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			     m_xkey.Doubles(),rans[1],mode);
}

void Resonance_Forward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

bool Resonance_Forward::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  m_spkey[2]=spkey[2];
  if (!m_zchannel) {
    m_spkey[0]=spkey[0];
    m_spkey[1]=spkey[3];
  }
  else {
    m_spkey[0]=spkey[3];
    m_spkey[1]=spkey[1];
    m_ykey[0]=ykey[0];
    m_ykey[1]=ykey[1];
    double logtau=.5*log(spkey[3]/spkey[2]);
    m_xkey[0]=logtau+ykey[2];
    m_xkey[2]=logtau-ykey[2];
    m_xkey[1]=0.;
    m_xkey[3]=0.;
  }
  return true;
}

Resonance_Backward::Resonance_Backward(const double mass,const double width,const double yexponent,
				     const std::string cinfo,Integration_Info *info): 
  m_mass(mass),
  m_width(width),
  m_yexponent(yexponent)
{
  char help[3];
  sprintf(help,"%i",int(mass*1000.));
  name=std::string(help);
  sprintf(help,"%i",int(width*1000.));
  name+=std::string(help);
  m_spkey.SetInfo(std::string("Resonance_")+std::string(help));
  name=std::string("Resonance_Backward");
  sprintf(help,"%i",int(100.*yexponent));
  m_ykey.SetInfo(std::string("Backward_")+std::string(help));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Resonance_Backward::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,int mode)
{
  CalculateLimits(spkey,ykey);
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]+=CE.DiceYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			      m_xkey.Doubles(),rans[1],mode);
}

void Resonance_Backward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_weight,1,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

bool Resonance_Backward::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  m_spkey[2]=spkey[2];
  if (!m_zchannel) {
    m_spkey[0]=spkey[0];
    m_spkey[1]=spkey[3];
  }
  else {
    m_spkey[0]=spkey[3];
    m_spkey[1]=spkey[1];
    m_ykey[0]=ykey[0];
    m_ykey[1]=ykey[1];
    double logtau=.5*log(spkey[3]/spkey[2]);
    m_xkey[0]=logtau+ykey[2];
    m_xkey[2]=logtau-ykey[2];
    m_xkey[1]=0.;
    m_xkey[3]=0.;
  }
  return true;
}

Resonance_Central::Resonance_Central(const double mass,const double width,
				   const std::string cinfo,Integration_Info *info): 
  m_mass(mass),
  m_width(width)
{
  char help[3];
  sprintf(help,"%i",int(mass*1000.));
  name=std::string(help);
  sprintf(help,"%i",int(width*1000.));
  name+=std::string(help);
  m_spkey.SetInfo(std::string("Resonance_")+std::string(help));
  name=std::string("Resonance_Central");
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Resonance_Central::GeneratePoint(Info_Key &spkey,Info_Key &ykey,const double *rans,int mode)
{
  CalculateLimits(spkey,ykey);
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]+=CE.DiceYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Resonance_Central::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

bool Resonance_Central::CalculateLimits(Info_Key &spkey,Info_Key &ykey) 
{
  m_spkey[2]=spkey[2];
  if (!m_zchannel) {
    m_spkey[0]=spkey[0];
    m_spkey[1]=spkey[3];
  }
  else {
    m_spkey[0]=spkey[3];
    m_spkey[1]=spkey[1];
    m_ykey[0]=ykey[0];
    m_ykey[1]=ykey[1];
    double logtau=.5*log(spkey[3]/spkey[2]);
    m_xkey[0]=logtau+ykey[2];
    m_xkey[2]=logtau-ykey[2];
    m_xkey[1]=0.;
    m_xkey[3]=0.;
  }
  return true;
}

