#include "Leading_Log_Z.H"

#include "Channel_Elements.H"
#include "Message.H"
#include "MyStrStream.H"

#include <stdio.h>

using namespace PHASIC;

Z_Channel_Base::Z_Channel_Base()
{
  rannum=2;
  p_vegas = new Vegas(2,250,name,0);
  rans = new double[2];
}

Z_Channel_Base::~Z_Channel_Base()
{
  delete [] rans;
  delete p_vegas;
}

double Z_Channel_Base::f(ATOOLS::Info_Key &key,double t)
{
  if (t<m_cut) return pow(m_cut,-m_beta);
  return pow(t,-m_beta);
}

double Z_Channel_Base::I(ATOOLS::Info_Key &key,double t)
{
  if (t<m_cut) return (t-key[0])*pow(m_cut,-m_beta);
  double I1=(m_cut-key[0])*pow(m_cut,-m_beta);
  return I1+(pow(t,1.-m_beta)-pow(m_cut,1.-m_beta))/(1.-m_beta);
}

double Z_Channel_Base::t(ATOOLS::Info_Key &key,double I)
{
  double I1=(m_cut-key[0])*pow(m_cut,-m_beta);
  if (I<I1) return I/pow(m_cut,-m_beta)+key[0];
  I-=I1;
  return pow(I*(1.-m_beta)+pow(m_cut,1.-m_beta),1.0/(1.-m_beta));
}

std::string Z_Channel_Base::ChID()
{
  return name;
}

void Z_Channel_Base::WriteOut(std::string pId) 
{ 
  p_vegas->WriteOut(pId); 
}

void Z_Channel_Base::ReadIn(std::string pId)   
{ 
  p_vegas->ReadIn(pId); 
}

void Z_Channel_Base::Optimize()  
{
  p_vegas->Optimize();
} 

void Z_Channel_Base::EndOptimize()  
{
  p_vegas->EndOptimize();
} 

Leading_Log_Z_QQ::Leading_Log_Z_QQ(const double beta,const double cut,
				   const std::string cinfo,
				   ATOOLS::Integration_Info *info)
{
  m_beta=beta;
  m_cut=cut;
  std::string help=ATOOLS::ToString(beta)+"_"+ATOOLS::ToString(cut);
  name=std::string("Leading_Log_Z_QQ_")+help;
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_Q_")+help);
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_Q_")+help);
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_QQ::
GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
	      const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short unsigned int i=0;i<2;++i)
    m_zkey[i][2]=1.0-t(m_zkey[i],I(m_zkey[i],m_zkey[i][0])*(1.0-rans[i])+
		       rans[i]*I(m_zkey[i],m_zkey[i][1]));
}

void Leading_Log_Z_QQ::GenerateWeight(const int mode) 
{
  for (short unsigned int i=0;i<2;++i)
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) 
      m_zkey[i]<<(I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]))/
	f(m_zkey[i],1.0-m_zkey[i][2]);
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_QQ::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  for (short unsigned int i=0;i<2;++i)
    rans[i]=(I(m_zkey[i],1.0-m_zkey[i][2])-I(m_zkey[i],m_zkey[i][0]))/
      (I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]));
  p_vegas->AddPoint(value,rans);
}

void Leading_Log_Z_QQ::CalculateLimits(ATOOLS::Info_Key &spkey,
				       ATOOLS::Info_Key &ykey) 
{
}

Leading_Log_Z_QG::
Leading_Log_Z_QG(const double beta,const double cut,
		 const std::string cinfo,ATOOLS::Integration_Info *info)
{
  m_beta=beta;
  m_cut=cut;
  std::string help=ATOOLS::ToString(beta)+"_"+ATOOLS::ToString(cut);
  name=std::string("Leading_Log_Z_QG_")+help;
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_Q_")+help);
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_G_")+help);
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_QG::
GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
	      const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short unsigned int i=0;i<2;++i) {
    double z=t(m_zkey[i],I(m_zkey[i],m_zkey[i][0])*(1.0-rans[i])+
	       rans[i]*I(m_zkey[i],m_zkey[i][1]));
    m_zkey[i][2]=i>0?z:(1.0-z);
  }
}

void Leading_Log_Z_QG::GenerateWeight(const int mode) 
{
  for (short unsigned int i=0;i<2;++i)
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double z=i>0?m_zkey[i][2]:(1.0-m_zkey[i][2]);
      m_zkey[i]<<(I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]))/
	f(m_zkey[i],z);
    }
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_QG::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  for (short unsigned int i=0;i<2;++i) {
    double z=i>0?m_zkey[i][2]:(1.0-m_zkey[i][2]);
    rans[i]=(I(m_zkey[i],z)-I(m_zkey[i],m_zkey[i][0]))/
      (I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]));
  }
  p_vegas->AddPoint(value,rans);
}

void Leading_Log_Z_QG::
CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}

Leading_Log_Z_GQ::
Leading_Log_Z_GQ(const double beta,const double cut,
		 const std::string cinfo,ATOOLS::Integration_Info *info)
{
  m_beta=beta;
  m_cut=cut;
  std::string help=ATOOLS::ToString(beta)+"_"+ATOOLS::ToString(cut);
  name=std::string("Leading_Log_Z_GQ_")+help;
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_G_")+help);
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_Q_")+help);
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_GQ::
GeneratePoint(ATOOLS::Info_Key &spkey,
	      ATOOLS::Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short unsigned int i=0;i<2;++i) {
    double z=t(m_zkey[i],I(m_zkey[i],m_zkey[i][0])*(1.0-rans[i])+
	       rans[i]*I(m_zkey[i],m_zkey[i][1]));
    m_zkey[i][2]=i<1?z:(1.0-z);
  }
}

void Leading_Log_Z_GQ::GenerateWeight(const int mode) 
{
  for (short unsigned int i=0;i<2;++i)
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double z=i<1?m_zkey[i][2]:(1.0-m_zkey[i][2]);
      m_zkey[i]<<(I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]))/
	f(m_zkey[i],z);
    }
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_GQ::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  for (short unsigned int i=0;i<2;++i) {
    double z=i<1?m_zkey[i][2]:(1.0-m_zkey[i][2]);
    rans[i]=(I(m_zkey[i],z)-I(m_zkey[i],m_zkey[i][0]))/
      (I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]));
  }
  p_vegas->AddPoint(value,rans);
}

void Leading_Log_Z_GQ::
CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}

Leading_Log_Z_GG::Leading_Log_Z_GG(const double beta,const double cut,
				   const std::string cinfo,
				   ATOOLS::Integration_Info *info)
{
  m_beta=beta;
  m_cut=cut;
  std::string help=ATOOLS::ToString(beta)+"_"+ATOOLS::ToString(cut);
  name=std::string("Leading_Log_Z_GG_")+help;
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_G_")+help);
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_G_")+help);
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_GG::
GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
	      const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short unsigned int i=0;i<2;++i) 
    m_zkey[i][2]=t(m_zkey[i],I(m_zkey[i],m_zkey[i][0])*(1.0-rans[i])+
		   rans[i]*I(m_zkey[i],m_zkey[i][1]));
}

void Leading_Log_Z_GG::GenerateWeight(const int mode) 
{
  for (short unsigned int i=0;i<2;++i)
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT)
      m_zkey[i]<<(I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]))/
	f(m_zkey[i],m_zkey[i][2]);
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_GG::AddPoint(double value)
{
  Single_Channel::AddPoint(value);
  for (short unsigned int i=0;i<2;++i) 
    rans[i]=(I(m_zkey[i],m_zkey[i][2])-I(m_zkey[i],m_zkey[i][0]))/
      (I(m_zkey[i],m_zkey[i][1])-I(m_zkey[i],m_zkey[i][0]));
  p_vegas->AddPoint(value,rans);
}

void Leading_Log_Z_GG::
CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}
