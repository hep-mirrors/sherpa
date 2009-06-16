#include "PHASIC++/Channels/Leading_Log.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <stdio.h>

using namespace PHASIC;

Leading_Log_Uniform::Leading_Log_Uniform(const double beta,const double factor,
					 const std::string cinfo,ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_factor(factor)
{
  name=std::string("Leading_Log_Uniform_")+ATOOLS::ToString((int)(100.*beta));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Leading_Log_Uniform::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rans,const int mode) 
{
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Uniform::GenerateWeight(const int mode) 
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
      double pole=m_spkey[2];
      if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

Leading_Log_Forward::Leading_Log_Forward(const double beta,const double factor,const double yexponent,
					 const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_beta(beta),
  m_factor(factor),
  m_yexponent(yexponent)
{
  name=std::string("Leading_Log_Forward_")+ATOOLS::ToString((int)(100.*beta));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta)+
		  std::string("_")+ATOOLS::ToString(yexponent));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Leading_Log_Forward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rans,const int mode) 
{
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Forward::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
      double pole=m_spkey[2];
      if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3]);
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

Leading_Log_Backward::Leading_Log_Backward(const double beta,const double factor,const double yexponent,
					   const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_beta(beta),
  m_factor(factor),
  m_yexponent(yexponent)
{
  name=std::string("Leading_Log_Backward_")+ATOOLS::ToString((int)(100.*beta));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta)+
		  std::string("_")+ATOOLS::ToString(yexponent));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Leading_Log_Backward::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					 const double *rans,int mode)
{
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,m_spkey[3]/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Backward::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
      double pole=m_spkey[2];
      if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,m_weight,m_spkey[0],m_spkey[1],m_spkey[3]);
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

Leading_Log_Central::Leading_Log_Central(const double beta,const double factor,
					 const std::string cinfo,ATOOLS::Integration_Info *info): 
  m_beta(beta),
  m_factor(factor)
{
  name=std::string("Leading_Log_Central_")+ATOOLS::ToString((int)(100.*beta));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
}

void Leading_Log_Central::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rans,int mode)
{
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Central::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
      double pole=m_spkey[2];
      if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3]);
    }
  }
  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral(m_spkey[3]/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),mode);
    }
  }
  weight=m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

