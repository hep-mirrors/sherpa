#include "PHASIC++/Channels/Resonance_Channels.H"
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

Resonance_RelicDensity::
Resonance_RelicDensity(const double mass,const double width,
		       const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width)
{
  m_name="Resonance_"+ATOOLS::ToString(mass)+"_RelicDensity";
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum  = 1;
  p_vegas = new Vegas(m_rannum,100,m_name,0);
  p_rans    = new double[m_rannum];
}

void Resonance_RelicDensity::GeneratePoint(const double *rns,const int mode)
{
  double * ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],p_rans[0]);
}

void Resonance_RelicDensity::GenerateWeight(const int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],
				       m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }
  p_rans[0] = m_sgridkey[0];
  double pw= p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()/m_spkey[2];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Resonance_Uniform::
Resonance_Uniform(const double mass,const double width,
		  const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width)
{
  m_name="Resonance_"+ATOOLS::ToString(mass)+"_Uniform";
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_rannum=2;
  p_vegas = new Vegas(2,100,m_name,0);
  p_rans  = new double[2];
}

void Resonance_Uniform::GeneratePoint(const double *rns,const int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYUniform(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),p_rans[1],mode);
}

void Resonance_Uniform::GenerateWeight(const int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],
				       m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYUniform(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),
				m_ygridkey[0],mode);
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


Resonance_Forward::
Resonance_Forward(const double mass,const double width,const double yexponent,
		  const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width),
  m_yexponent(yexponent)
{
  m_name="Resonance_"+ATOOLS::ToString(mass)+"_Forward_"+ATOOLS::ToString(yexponent);
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
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
  p_vegas = new Vegas(2,100,m_name,0);
  p_rans  = new double[2];
}

void Resonance_Forward::GeneratePoint(const double *rns,const int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYForward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),p_rans[1],mode);
}

void Resonance_Forward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],
				       m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYForward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
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


Resonance_Backward::
Resonance_Backward(const double mass,const double width,const double yexponent,
		   const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width),
  m_yexponent(yexponent)
{
  m_name="Resonance_"+ATOOLS::ToString(mass)+"_Backward_"+ATOOLS::ToString(yexponent);
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
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
  p_vegas = new Vegas(2,100,m_name,0);
  p_rans  = new double[2];
}

void Resonance_Backward::GeneratePoint(const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYBackward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),p_rans[1],mode);
}

void Resonance_Backward::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],
				       m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYBackward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],mode);
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


Resonance_Central::
Resonance_Central(const double mass,const double width,
		  const std::string cinfo,ATOOLS::Integration_Info *info,int mode):
  ISR_Channel_Base(info), m_mass(mass), m_width(width)
{
  m_name="Resonance_"+ATOOLS::ToString(mass)+"_Central";
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_rannum=1;
  if (mode==3) m_rannum=2;
  p_vegas = new Vegas(m_rannum,100,m_name,0);
  p_rans  = new double[2];
}

void Resonance_Central::GeneratePoint(const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  p_rans[0]=ran[0];
  if (mode==3) p_rans[1]=ran[1];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYCentral(sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),p_rans[1],mode);
}

void Resonance_Central::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],
				       m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
      m_ykey<<CE.WeightYCentral(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),
				m_ygridkey[0],mode);
    }
  }
  p_rans[0] = m_sgridkey[0];
  p_rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}



