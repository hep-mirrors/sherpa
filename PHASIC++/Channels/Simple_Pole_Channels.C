#include "PHASIC++/Channels/Simple_Pole_Channels.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Simple_Pole_RelicDensity::
Simple_Pole_RelicDensity(const double exponent, const std::string cinfo,
				ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_exponent(exponent)
{
  m_name="Simple_Pole_"+ATOOLS::ToString(exponent)+"_RelicDensity";
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_zchannel = m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum   = 1;
  p_vegas    = new Vegas(m_rannum,100,m_name,0);
  p_rans     = new double[m_rannum];
	msg_Out()<<"exponent="<<m_exponent<<"\n"; //debugging
}

void Simple_Pole_RelicDensity::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
	m_spkey[3] = CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],p_rans[0]);
}

void Simple_Pole_RelicDensity::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
			m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3],
							m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  p_rans[0] = m_sgridkey[0];
  double pw = p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight();
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Simple_Pole_DM_Annihilation::
Simple_Pole_DM_Annihilation(const double exponent, const double mass1, const double mass2,
			    const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),m_exponent(exponent)
{
  m_mass[0] = mass1;
  m_mass[1] = mass2;
  m_name="Simple_Pole_"+ATOOLS::ToString(exponent)+"_DM_Annihilation";
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),3,0,info);
  m_cosxikey.Assign(cinfo+std::string("::cosXi"),3,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_xgridkey.Assign(m_xkey.Info(),1,0,info);
  m_cosgridkey.Assign(m_cosxikey.Info(),1,0,info);
  m_zchannel = m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum   = 2;
  p_vegas    = new Vegas(m_rannum,100,m_name,0);
  p_rans     = new double[m_rannum];
}

void Simple_Pole_DM_Annihilation::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<m_rannum;i++) p_rans[i]=ran[i];
  m_spkey[3] = CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],p_rans[0]);

  // for now, all p_rans[0]. Change to [1] and [2] when m_rannum fixed
	m_cosxikey[2] = CE.GenerateDMAngleUniform(p_rans[1],3);
  m_xkey[2] = CE.GenerateDMRapidityUniform(m_mass,m_spkey.Doubles(),m_xkey.Doubles(),
					   m_cosxikey[2], p_rans[0], 3);
}

void Simple_Pole_DM_Annihilation::GenerateWeight(const int & mode)
{
  // this needs looking at
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3],
					m_sgridkey[0]);
    }
  }
  if (m_spkey[4]>0.0) { p_vegas->ConstChannel(0); m_spkey<<M_PI*2.0; }

  p_rans[0] = m_sgridkey[0];
  double pw = p_vegas->GenerateWeight(p_rans);
  m_weight  = pw*m_spkey.Weight()/m_spkey[3];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

Simple_Pole_Uniform::Simple_Pole_Uniform(const double exponent,const std::string cinfo,
					 ATOOLS::Integration_Info *info,
					 const size_t mode):
  ISR_Channel_Base(info),
  m_exponent(exponent),
  m_mode(mode)
{
  m_name="Simple_Pole_"+ATOOLS::ToString(exponent)+"_Uniform";
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum = 2;
  p_vegas  = new Vegas(2,100,m_name,0);
  p_rans   = new double[2];
}

void Simple_Pole_Uniform::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYUniform(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),
				p_rans[1],m_mode);
}

void Simple_Pole_Uniform::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3],
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

Simple_Pole_Forward::Simple_Pole_Forward(const double sexponent,const double yexponent,
					 const std::string cinfo,ATOOLS::Integration_Info *info,
					 const size_t mode):
  ISR_Channel_Base(info),
  m_sexponent(sexponent),
  m_yexponent(yexponent),
  m_mode(mode)
{
  m_name="Simple_Pole_"+ATOOLS::ToString(sexponent)+"_Forward_"+ATOOLS::ToString(yexponent);
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(sexponent));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
	m_kp1key.Assign("k_perp_1",4,1,info);
	m_kp2key.Assign("k_perp_2",4,1,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum=2;
  p_vegas = new Vegas(2,100,m_name,0);
  p_rans  = new double[2];
}

void Simple_Pole_Forward::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_sexponent,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYForward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),p_rans[1],m_mode);
}

void Simple_Pole_Forward::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_sexponent,m_spkey[0],m_spkey[1],m_spkey[3],
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

Simple_Pole_Backward::Simple_Pole_Backward(const double sexponent,const double yexponent,
					   const std::string cinfo,ATOOLS::Integration_Info *info,
					   const size_t mode) :
  ISR_Channel_Base(info),
  m_sexponent(sexponent),
  m_yexponent(yexponent),
  m_mode(mode)
{
  m_name="Simple_Pole_"+ATOOLS::ToString(sexponent)+"_Backward_"+ATOOLS::ToString(yexponent);
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(sexponent));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum=2;
  p_vegas = new Vegas(2,100,m_name,0);
  p_rans  = new double[2];
}

void Simple_Pole_Backward::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_sexponent,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYBackward(m_yexponent,sred/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),p_rans[1],m_mode);
}

void Simple_Pole_Backward::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_sexponent,m_spkey[0],m_spkey[1],m_spkey[3],
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

Simple_Pole_Central::Simple_Pole_Central(const double exponent,const std::string cinfo,
					 ATOOLS::Integration_Info *info,
					 const size_t mode):
  ISR_Channel_Base(info),
  m_exponent(exponent),
  m_mode(mode)
{
  m_name="Simple_Pole_"+ATOOLS::ToString(exponent)+"_Central";
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(cinfo+std::string("::s'"),5,0,info);
  m_ykey.Assign(cinfo+std::string("::y"),3,0,info);
  m_xkey.Assign(cinfo+std::string("::x"),5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  m_rannum=2;
  p_vegas = new Vegas(m_rannum,100,m_name,0);
  p_rans  = new double[2];
}

void Simple_Pole_Central::GeneratePoint(const double *rns)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) p_rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],p_rans[0]);
  double sred = SelectS(m_spkey[3],m_spkey[4])-(m_kp1key(0)+m_kp2key(0)).Abs2();
  m_ykey[2]=CE.GenerateYCentral(sred/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),p_rans[1],
				m_mode);
}

void Simple_Pole_Central::GenerateWeight(const int & mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3],
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
  double pw = p_vegas->GenerateWeight(p_rans);
  m_weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}
