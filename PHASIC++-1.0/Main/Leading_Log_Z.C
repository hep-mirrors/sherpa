#include "Leading_Log_Z.H"

#include "Channel_Elements.H"
#include "Message.H"

#include <stdio.h>

using namespace PHASIC;

Leading_Log_Z_QQ::Leading_Log_Z_QQ(const double beta,const double cut,
				   const std::string cinfo,ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_cut(cut)
{
  char help[3];
  sprintf(help,"%i_%i",int(beta*100.),int(cut*100.));
  name=std::string("Leading_Log_Z_QQ_")+std::string(help);
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_Q_")+std::string(help));
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_Q_")+std::string(help));
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_QQ::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short int i=0;i<2;++i) {
    double zc=m_zkey[i][1]*m_cut, rc=(zc-m_zkey[i][0])/(m_zkey[i][1]-m_zkey[i][0]);
    if (rans[i]<rc) {
      m_zkey[i][2]=CE.LLPropMomenta(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,rans[i]/rc);
    }
    else {
      m_zkey[i][2]=m_zkey[i][1]*(m_cut+(1.-m_cut)*(1.-rans[i])/(1.-rc));
    }
  }
}

void Leading_Log_Z_QQ::GenerateWeight(const int mode) 
{
  weight=0.;
  for (short int i=0;i<2;++i) {
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      if (m_zkey[i][2]>=m_zkey[i][0] && m_zkey[i][2]<=m_zkey[i][1]) {
	double zc=m_zkey[i][1]*m_cut;
	if (m_zkey[i][2]<zc) {
	  m_zkey[i]<<1./CE.LLPropWeight(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,m_zkey[i][2]);
	}
	else {
	  m_zkey[i]<<m_zkey[i][1]-zc;
	}
      }
    }
  }
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_QQ::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}

Leading_Log_Z_QG::Leading_Log_Z_QG(const double beta,const double cut,
				   const std::string cinfo,ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_cut(cut)
{
  char help[3];
  sprintf(help,"%i_%i",int(beta*100.),int(cut*100.));
  name=std::string("Leading_Log_Z_QG_")+std::string(help);
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_Q_")+std::string(help));
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_G_")+std::string(help));
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_QG::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short int i=0;i<2;++i) {
    double zc=m_zkey[i][1]*m_cut, rc=(zc-m_zkey[i][0])/(m_zkey[i][1]-m_zkey[i][0]);
    if (rans[i]<rc) {
      double value=CE.LLPropMomenta(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,rans[i]/rc);
      if (i==0) m_zkey[i][2]=value;
      else m_zkey[i][2]=1.-value;
    }
    else {
      double value=m_zkey[i][1]*(m_cut+(1.-m_cut)*(1.-rans[i])/(1.-rc));
      if (i==0) m_zkey[i][2]=value;
      else m_zkey[i][2]=1.-value;
    }
  }
}

void Leading_Log_Z_QG::GenerateWeight(const int mode) 
{
  weight=0.;
  for (short int i=0;i<2;++i) {
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double value=m_zkey[i][2];
      if (i==1) value=1.-value;
      if (value>=m_zkey[i][0] && value<=m_zkey[i][1]) {
	double zc=m_zkey[i][1]*m_cut;
	if (value<zc) {
	  m_zkey[i]<<1./CE.LLPropWeight(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,value);
	}
	else {
	  m_zkey[i]<<m_zkey[i][1]-zc;
	}
      }
    }
  }
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_QG::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}

Leading_Log_Z_GQ::Leading_Log_Z_GQ(const double beta,const double cut,
				   const std::string cinfo,ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_cut(cut)
{
  char help[3];
  sprintf(help,"%i_%i",int(beta*100.),int(cut*100.));
  name=std::string("Leading_Log_Z_GQ_")+std::string(help);
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_G_")+std::string(help));
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_Q_")+std::string(help));
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_GQ::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short int i=0;i<2;++i) {
    double zc=m_zkey[i][1]*m_cut, rc=(zc-m_zkey[i][0])/(m_zkey[i][1]-m_zkey[i][0]);
    if (rans[i]<rc) {
      double value=CE.LLPropMomenta(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,rans[i]/rc);
      if (i==1) m_zkey[i][2]=value;
      else m_zkey[i][2]=1.-value;
    }
    else {
      double value=m_zkey[i][1]*(m_cut+(1.-m_cut)*(1.-rans[i])/(1.-rc));
      if (i==1) m_zkey[i][2]=value;
      else m_zkey[i][2]=1.-value;
    }
  }
}

void Leading_Log_Z_GQ::GenerateWeight(const int mode) 
{
  weight=0.;
  for (short int i=0;i<2;++i) {
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double value=m_zkey[i][2];
      if (i==0) value=1.-value;
      if (value>=m_zkey[i][0] && value<=m_zkey[i][1]) {
	double zc=m_zkey[i][1]*m_cut;
	if (value<zc) {
	  m_zkey[i]<<1./CE.LLPropWeight(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,value);
	}
	else {
	  m_zkey[i]<<m_zkey[i][1]-zc;
	}
      }
    }
  }
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_GQ::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}

Leading_Log_Z_GG::Leading_Log_Z_GG(const double beta,const double cut,
				   const std::string cinfo,ATOOLS::Integration_Info *info):
  m_beta(beta),
  m_cut(cut)
{
  char help[3];
  sprintf(help,"%i_%i",int(beta*100.),int(cut*100.));
  name=std::string("Leading_Log_Z_GG_")+std::string(help);
  m_zkey[0].SetInfo(std::string("Leading_Log_Z_G_")+std::string(help));
  m_zkey[1].SetInfo(std::string("Leading_Log_Z_G_")+std::string(help));
  m_zkey[0].Assign(std::string("z_1")+cinfo,3,0,info);
  m_zkey[1].Assign(std::string("z_2")+cinfo,3,0,info);
}

void Leading_Log_Z_GG::GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode) 
{
  CalculateLimits(spkey,ykey);
  for (short int i=0;i<2;++i) {
    double zc=m_zkey[i][1]*m_cut, rc=(zc-m_zkey[i][0])/(m_zkey[i][1]-m_zkey[i][0]);
    if (rans[i]<rc) {
      m_zkey[i][2]=1.-CE.LLPropMomenta(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,rans[i]/rc);
    }
    else {
      m_zkey[i][2]=1.-m_zkey[i][1]*(m_cut+(1.-m_cut)*(1.-rans[i])/(1.-rc));
    }
  }
}

void Leading_Log_Z_GG::GenerateWeight(const int mode) 
{
  weight=0.;
  for (short int i=0;i<2;++i) {
    if (m_zkey[i].Weight()==ATOOLS::UNDEFINED_WEIGHT) {
      double value=1.-m_zkey[i][2];
      if (value>=m_zkey[i][0] && value<=m_zkey[i][1]) {
	double zc=m_zkey[i][1]*m_cut;
	if (value<zc) {
	  m_zkey[i]<<1./CE.LLPropWeight(1.-m_beta,m_zkey[i][1],m_zkey[i][0],zc,value);
	}
	else {
	  m_zkey[i]<<m_zkey[i][1]-zc;
	}
      }
    }
  }
  weight=m_zkey[0].Weight()*m_zkey[1].Weight();
}

void Leading_Log_Z_GG::CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey) 
{
}
