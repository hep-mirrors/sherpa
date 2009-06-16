#include "PHASIC++/Channels/ISR_Vegas.H"

#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <stdio.h>

using namespace PHASIC;

ISR_Channel_Base::ISR_Channel_Base(ATOOLS::Integration_Info *info)
{
  m_kp1key.Assign("k_perp_1",4,1,info);
  m_kp2key.Assign("k_perp_2",4,1,info);
}

ISR_Channel_Base::~ISR_Channel_Base()
{
  delete p_vegas;
}

std::string ISR_Channel_Base::ChID()
{
  return name;
}

void ISR_Channel_Base::WriteOut(std::string pId) 
{ 
  p_vegas->WriteOut(pId); 
}

void ISR_Channel_Base::ReadIn(std::string pId)   
{ 
  p_vegas->ReadIn(pId); 
}

void ISR_Channel_Base::Optimize()  
{
  p_vegas->Optimize();
} 

void ISR_Channel_Base::EndOptimize()  
{
  p_vegas->EndOptimize();
} 

Threshold_Uniform_V::Threshold_Uniform_V(const double mass,const double sexp,const std::string cinfo,
				     ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_mass(mass), m_sexp(sexp)
{
  name=std::string("Threshold_Uniform_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Threshold_Uniform_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.ThresholdMomenta(m_sexp,m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Uniform_V::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_sexp,m_mass,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Threshold_Uniform_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}

Threshold_Forward_V::Threshold_Forward_V(const double mass,const double sexp,const double yexponent,
				     const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_mass(mass), m_sexp(sexp), 
  m_yexponent(yexponent)
{
  name=std::string("Threshold_Forward_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Threshold_Forward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.ThresholdMomenta(m_sexp,m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Forward_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_sexp,m_mass,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Threshold_Forward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}

Threshold_Backward_V::Threshold_Backward_V(const double mass,const double sexp,const double yexponent,
				       const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_mass(mass), m_sexp(sexp), 
  m_yexponent(yexponent)
{
  name=std::string("Threshold_Backward_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Threshold_Backward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				       const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.ThresholdMomenta(m_sexp,m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Backward_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_sexp,m_mass,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],mode);
     }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Threshold_Backward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}

Threshold_Central_V::Threshold_Central_V(const double mass,const double sexp,const std::string cinfo,
				     ATOOLS::Integration_Info *info,int mode):
  ISR_Channel_Base(info),
  m_mass(mass), m_sexp(sexp)
{
  name=std::string("Threshold_Central_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Threshold_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=1;
  if (mode==3) rannum=2;
  p_vegas = new Vegas(rannum,100,name,0);
  rans  = new double[2];
}

void Threshold_Central_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  rans[0]=ran[0];
  if (mode==3) rans[1]=ran[1];
  m_spkey[3]=CE.ThresholdMomenta(m_sexp,m_mass,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Threshold_Central_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.ThresholdWeight(m_sexp,m_mass,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Threshold_Central_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Resonance_Uniform_V::Resonance_Uniform_V(const double mass,const double width,
				     const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width)
{
  name=std::string("Resonance_Uniform_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Resonance_Uniform_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Resonance_Uniform_V::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Resonance_Uniform_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Resonance_Forward_V::Resonance_Forward_V(const double mass,const double width,const double yexponent,
				   const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width),
  m_yexponent(yexponent)
{
  name=std::string("Resonance_Forward_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Resonance_Forward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Resonance_Forward_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Resonance_Forward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Resonance_Backward_V::Resonance_Backward_V(const double mass,const double width,const double yexponent,
				       const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width),
  m_yexponent(yexponent)
{
  name=std::string("Resonance_Backward_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Resonance_Backward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				       const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Resonance_Backward_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Resonance_Backward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Resonance_Central_V::Resonance_Central_V(const double mass,const double width,
				     const std::string cinfo,ATOOLS::Integration_Info *info,int mode): 
  ISR_Channel_Base(info),
  m_mass(mass),
  m_width(width)
{
  name=std::string("Resonance_Central_")+ATOOLS::ToString((int)(100.*mass+0.01));
  m_spkey.SetInfo(std::string("Resonance_")+ATOOLS::ToString(mass));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=1;
  if (mode==3) rannum=2;
  p_vegas = new Vegas(rannum,100,name,0);
  rans  = new double[2];
}

void Resonance_Central_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
				      const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  rans[0]=ran[0];
  if (mode==3) rans[1]=ran[1];
  m_spkey[3]=CE.MassivePropMomenta(m_mass,m_width,1,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Resonance_Central_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MassivePropWeight(m_mass,m_width,1,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Resonance_Central_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Simple_Pole_Uniform_V::Simple_Pole_Uniform_V(const double exponent,const std::string cinfo,
					 ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_exponent(exponent)
{
  name=std::string("Simple_Pole_Uniform_")+ATOOLS::ToString((int)(100.*exponent+0.01));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Simple_Pole_Uniform_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Uniform_V::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Simple_Pole_Uniform_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Simple_Pole_Forward_V::Simple_Pole_Forward_V(const double sexponent,const double yexponent,
					 const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_sexponent(sexponent), 
  m_yexponent(yexponent)
{
  name=std::string("Simple_Pole_Forward_")+ATOOLS::ToString((int)(100.*sexponent+0.01));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(sexponent));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Simple_Pole_Forward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_sexponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Forward_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_sexponent,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Simple_Pole_Forward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Simple_Pole_Backward_V::Simple_Pole_Backward_V(const double sexponent,const double yexponent,
					   const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_sexponent(sexponent), 
  m_yexponent(yexponent)
{
  name=std::string("Simple_Pole_Backward_")+ATOOLS::ToString((int)(100.*sexponent+0.01));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(sexponent));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Simple_Pole_Backward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					 const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.MasslessPropMomenta(m_sexponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Backward_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_sexponent,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Simple_Pole_Backward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Simple_Pole_Central_V::Simple_Pole_Central_V(const double exponent,const std::string cinfo,
					 ATOOLS::Integration_Info *info,int mode):
  ISR_Channel_Base(info),
  m_exponent(exponent)
{
  name=std::string("Simple_Pole_Central_")+ATOOLS::ToString((int)(100.*exponent+0.01));
  m_spkey.SetInfo(std::string("Simple_Pole_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=1;
  if (mode==3) rannum=2;
  p_vegas = new Vegas(rannum,100,name,0);
  rans  = new double[2];
}

void Simple_Pole_Central_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  rans[0]=ran[0];
  if (mode==3) rans[1]=ran[1];
  m_spkey[3]=CE.MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Simple_Pole_Central_V::GenerateWeight(int mode)
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Simple_Pole_Central_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Leading_Log_Uniform_V::Leading_Log_Uniform_V(const double beta,const double factor,
					 const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_beta(beta),
  m_factor(factor)
{
  name=std::string("Leading_Log_Uniform_")+ATOOLS::ToString((int)(100.*beta+0.01));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Leading_Log_Uniform_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Uniform_V::GenerateWeight(const int mode) 
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
    double pole=m_spkey[2];
    if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Leading_Log_Uniform_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Leading_Log_Forward_V::Leading_Log_Forward_V(const double beta,const double factor,const double yexponent,
					 const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_beta(beta),
  m_factor(factor),
  m_yexponent(yexponent)
{
  name=std::string("Leading_Log_Forward_")+ATOOLS::ToString((int)(100.*beta+0.01));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta));
  m_ykey.SetInfo(std::string("Forward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Leading_Log_Forward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Forward_V::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
    double pole=m_spkey[2];
    if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Leading_Log_Forward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Leading_Log_Backward_V::Leading_Log_Backward_V(const double beta,const double factor,const double yexponent,
					   const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
  m_beta(beta),
  m_factor(factor),
  m_yexponent(yexponent)
{
  name=std::string("Leading_Log_Backward_")+ATOOLS::ToString((int)(100.*beta+0.01));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta));
  m_ykey.SetInfo(std::string("Backward_")+ATOOLS::ToString(yexponent));
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Leading_Log_Backward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					 const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Backward_V::GenerateWeight(int mode)
{
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
    double pole=m_spkey[2];
    if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void Leading_Log_Backward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Leading_Log_Central_V::Leading_Log_Central_V(const double beta,const double factor,
					 const std::string cinfo,ATOOLS::Integration_Info *info,int mode): 
  ISR_Channel_Base(info),
  m_beta(beta),
  m_factor(factor)
{
  name=std::string("Leading_Log_Central_")+ATOOLS::ToString((int)(100.*beta+0.01));
  m_spkey.SetInfo(std::string("Leading_Log_")+ATOOLS::ToString(beta));
  m_ykey.SetInfo("Central");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=1;
  if (mode==3) rannum=2;
  p_vegas = new Vegas(rannum,100,name,0);
  rans  = new double[2];
}

void Leading_Log_Central_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  rans[0]=ran[0];
  if (mode==3) rans[1]=ran[1];
  double pole=m_spkey[2];
  if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
  m_spkey[3]=CE.LLPropMomenta(1.-m_beta,pole,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Leading_Log_Central_V::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    if (m_spkey[3]<m_spkey[0] || m_spkey[3]>m_spkey[1]) return;
    double pole=m_spkey[2];
    if (ATOOLS::IsEqual(m_spkey[2],m_spkey[1])) pole*=m_factor;
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(1.-m_beta,pole,m_spkey[0],m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Leading_Log_Central_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


LBS_Compton_Peak_Uniform_V::LBS_Compton_Peak_Uniform_V(const double exponent,const double pole,
						   const std::string cinfo,ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
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
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void LBS_Compton_Peak_Uniform_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					     const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Uniform_V::GenerateWeight(const int mode) 
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void LBS_Compton_Peak_Uniform_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


LBS_Compton_Peak_Forward_V::LBS_Compton_Peak_Forward_V(const double exponent,const double pole,
						   const double yexponent,
						   const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
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
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void LBS_Compton_Peak_Forward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					     const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			     m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Forward_V::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYForward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void LBS_Compton_Peak_Forward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


LBS_Compton_Peak_Backward_V::LBS_Compton_Peak_Backward_V(const double exponent,const double pole,
						     const double yexponent,
						     const std::string cinfo,ATOOLS::Integration_Info *info): 
  ISR_Channel_Base(info),
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
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void LBS_Compton_Peak_Backward_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					      const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
			      m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Backward_V::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYBackward(m_yexponent,(m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				 m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
} 

void LBS_Compton_Peak_Backward_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


LBS_Compton_Peak_Central_V::LBS_Compton_Peak_Central_V(const double exponent,const double pole,
						   const std::string cinfo,ATOOLS::Integration_Info *info,int mode): 
  ISR_Channel_Base(info),
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
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=1;
  if (mode==3) rannum=2;
  p_vegas = new Vegas(rannum,100,name,0);
  rans  = new double[2];
}

void LBS_Compton_Peak_Central_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
					     const double *rns,int mode)
{
  double *ran = p_vegas->GeneratePoint(rns);
  rans[0]=ran[0];
  if (mode==3) rans[1]=ran[1];
  double help=CE.LLPropMomenta(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],rans[0]);
  if (m_spkey[0]<m_spkey[2]*m_pole && m_spkey[2]*m_pole<m_spkey[1]) {
    m_spkey[3]=help-m_spkey[1]+m_spkey[2]*m_pole;
    if (m_spkey[3]<m_spkey[0]) m_spkey[3]=help+(m_spkey[2]*m_pole-m_spkey[0]);
  }
  else {
    m_spkey[3]=help;
  }
  m_ykey[2]=CE.DiceYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void LBS_Compton_Peak_Central_V::GenerateWeight(int mode)
{
  weight=0.;
  if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
    double help=m_spkey[3];
    if (m_spkey[0]<m_spkey[2]*m_pole || m_spkey[2]*m_pole<m_spkey[1]) {
      if (m_spkey[3]>m_pole*m_spkey[2]) help=m_spkey[3]-(m_spkey[2]*m_pole-m_spkey[0]);
      else help=m_spkey[3]+m_spkey[1]-m_spkey[2]*m_pole;
    }
    if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      m_spkey<<1./CE.LLPropWeight(m_exponent,m_spkey[2],m_spkey[0],m_spkey[1],help,m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYCentral((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void LBS_Compton_Peak_Central_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}


Flat_ISR_V::Flat_ISR_V(const double exponent,const std::string cinfo,
		       ATOOLS::Integration_Info *info):
  ISR_Channel_Base(info),
  m_exponent(exponent)
{
  name=std::string("Simple_Pole_PDF_Uniform_")+ATOOLS::ToString((int)(100.*exponent+0.01));
  m_spkey.SetInfo(std::string("Simple_Pole_PDF_")+ATOOLS::ToString(exponent));
  m_ykey.SetInfo("Uniform");
  m_spkey.Assign(std::string("s'")+cinfo,4,0,info);
  m_ykey.Assign(std::string("y")+cinfo,3,0,info);
  m_xkey.Assign(std::string("x")+cinfo,5,0,info);
  m_sgridkey.Assign(m_spkey.Info(),1,0,info);
  m_ygridkey.Assign(m_ykey.Info(),1,0,info);
  m_zchannel=m_spkey.Name().find("z-channel")!=std::string::npos;
  rannum=2;
  p_vegas = new Vegas(2,100,name,0);
  rans  = new double[2];
}

void Flat_ISR_V::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
			       const double *rns,const int mode) 
{
  double *ran = p_vegas->GeneratePoint(rns);
  for(int i=0;i<2;i++) rans[i]=ran[i];
  m_spkey[3]=CE.FlatMomenta(m_exponent,m_spkey[0],m_spkey[1],rans[0]);
  m_ykey[2]=CE.DiceYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],
			    m_xkey.Doubles(),m_ykey.Doubles(),rans[1],mode);
}

void Flat_ISR_V::GenerateWeight(const int mode) 
{
  if (m_spkey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_spkey[3]>=m_spkey[0] && m_spkey[3]<=m_spkey[1]) {
      m_spkey<<1./CE.FlatWeight(m_exponent,m_spkey[0],
				m_spkey[1],m_spkey[3],m_sgridkey[0]);
    }
  }

  if (m_ykey.Weight()==ATOOLS::UNDEFINED_WEIGHT) {
    if (m_ykey[2]>=m_ykey[0] && m_ykey[2]<=m_ykey[1]) {
      m_ykey<<CE.WeightYUniform((m_spkey[3]-(m_kp1key(0)+m_kp2key(0)).Abs2())/m_spkey[2],m_xkey.Doubles(),
				m_ykey.Doubles(),m_ygridkey[0],mode);
    }
  }
  rans[0] = m_sgridkey[0];
  rans[1] = m_ygridkey[0];
  double pw= p_vegas->GenerateWeight(rans);
  weight=pw*m_spkey.Weight()*m_ykey.Weight()/m_spkey[2];
}

void Flat_ISR_V::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  p_vegas->AddPoint(value,rans);
}

