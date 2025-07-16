#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


S1Channel::S1Channel(int _nin,int _nout,Flavour * fl,Flavour res) :
  Single_Channel(_nin,_nout,fl)
{  
  if (m_nin != 2 || m_nout!=2) {
    msg_Error()<<"Tried to initialize S1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    Abort();
  }
  m_rannum = 2;
  delete p_rans;
  p_rans   = new double[m_rannum];
  m_name   = "S-Channel";

  s      = smax  = pt2max = sqr(ATOOLS::rpa->gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);

  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(m_rannum,100,m_name);
}

S1Channel::~S1Channel()
{
  delete p_vegas;
}

void S1Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran=0) {
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=1.0;
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+p_ms[2]-p_ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.Isotropic2Momenta(p[0]+p[1],p_ms[2],p_ms[3],p[2],p[3],ran[0],ran[1],-ctmax,ctmax);
}

void S1Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) {
  double ctmax=1.0;
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+p_ms[2]-p_ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  double rans[2];
  m_weight = 1. / ( CE.Isotropic2Weight(p[2],p[3],rans[0],rans[1],-ctmax,ctmax) *
		    pow(2.*M_PI,2.*3.-4.) );
  m_weight *= p_vegas->GenerateWeight(rans);
}

void S1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

std::string S1Channel::ChID() 
{
  return std::string("S-Channel");
}

namespace PHASIC {

  class S1_Channel_Generator: public Channel_Generator {
  public:
    
    S1_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new S1Channel(p_proc->NIn(),p_proc->NOut(),
				  (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }

  };// end of class S1_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(S1_Channel_Generator,"SChannel",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,S1_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new S1_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    S1_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"2->2 S-channel integrator";
}

T1Channel::T1Channel(int _nin,int _nout,Flavour * fl,Flavour res) :
    Single_Channel(_nin,_nout,fl)
{  
  if (m_nout != 2 || m_nin!=2) {
    msg_Error()<<"Tried to initialize T1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    Abort();
  }
  m_rannum = 2;
  delete p_rans;
  p_rans   = new double[m_rannum];
  s      = smax  = pt2max = sqr(ATOOLS::rpa->gen.Ecms());
  pt2min = 0.0;
  E      = 0.5 * sqrt(s);
  m_name = "T-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(m_rannum,100,m_name);
}

T1Channel::~T1Channel()
{
  delete p_vegas;
}

void T1Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran =0) 
{
  double ctmax=1.0;
  double *ran = p_vegas->GeneratePoint(_ran);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+p_ms[2]-p_ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.TChannelMomenta(p[0],p[1],p[2],p[3],p_ms[2],p_ms[3],0.,
		     .5,ctmax,-ctmax,ran[0],ran[1]);
}

void T1Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) 
{
  double ctmax=1.0;
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+p_ms[2]-p_ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  double rans[2];
  m_weight = 1. / ( CE.TChannelWeight(p[0],p[1],p[2],p[3],0.,
				    .5,ctmax,-ctmax,rans[0],rans[1]) 
		  * pow(2.*M_PI,2*3.-4.) );
  m_weight *= p_vegas->GenerateWeight(rans);
}

void T1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}

namespace PHASIC {

  class T1_Channel_Generator: public Channel_Generator {
  public:
    
    T1_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new T1Channel(p_proc->NIn(),p_proc->NOut(),
				  (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }

  };// end of class T1_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(T1_Channel_Generator,"TChannel",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,T1_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new T1_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    T1_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"2->2 T-channel integrator";
}

U1Channel::U1Channel(int _nin,int _nout,Flavour * fl,Flavour res) :
    Single_Channel(_nin,_nout,fl)
{  
  if (m_nout != 2 || m_nin!=2) {
    msg_Error()<<"Tried to initialize U1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    Abort();
  }
  m_rannum = 2;
  delete p_rans;
  p_rans   = new double[m_rannum];
  s      = smax  = pt2max = sqr(ATOOLS::rpa->gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  m_name = "U-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(m_rannum,100,m_name);
}

U1Channel::~U1Channel()
{
  delete p_vegas;
}

void U1Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran =0) 
{
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=1.0;
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+p_ms[2]-p_ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.TChannelMomenta(p[0],p[1],p[3],p[2],p_ms[3],p_ms[2],0.,
		     0.5,ctmax,-ctmax,ran[0],ran[1]);
}

void U1Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) 
{
  double ctmax=1.0;
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+p_ms[2]-p_ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  double rans[2];
  m_weight = 1. / ( CE.TChannelWeight(p[0],p[1],p[3],p[2],0.,
				    .5,ctmax,-ctmax,rans[0],rans[1]) 
		  * pow(2.*M_PI,2*3.-4.) );
  m_weight *= p_vegas->GenerateWeight(rans);
}

void U1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}

std::string U1Channel::ChID() 
{
  return std::string("U-Channel");
}

namespace PHASIC {

  class U1_Channel_Generator: public Channel_Generator {
  public:
    
    U1_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new U1Channel(p_proc->NIn(),p_proc->NOut(),
				  (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }

  };// end of class U1_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(U1_Channel_Generator,"UChannel",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,U1_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new U1_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    U1_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"2->2 U-channel integrator";
}

Decay2Channel::Decay2Channel(int _nin,int _nout,const Flavour * fl,Flavour res) 
{
  m_nin = _nin; m_nout = _nout;
  p_ms[0] = ATOOLS::sqr(res.Mass());
  for (int i=m_nin;i<m_nin+m_nout;i++) {
    p_ms[i] = ATOOLS::sqr(fl[i-1].Mass());
  }
  if (m_nout != 2 || m_nin!=1) {
    msg_Error()<<"Tried to initialize Decay2Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    Abort();
  }
  m_rannum = 2;
  delete p_rans;
  p_rans   = new double[m_rannum];
  s      = smax  = pt2max = sqr(ATOOLS::rpa->gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  m_name = "Decay2-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void Decay2Channel::GeneratePoint(ATOOLS::Vec4D * p,double * _ran=0) {
  CE.Isotropic2Momenta(p[0],p_ms[1],p_ms[2],p[1],p[2],_ran[0],_ran[1],-1.,1.);
}

void Decay2Channel::GenerateWeight(ATOOLS::Vec4D * p) {
  double d1, d2;
  m_weight = 1. / ( CE.Isotropic2Weight(p[1],p[2],d1,d2,-1.,1.) * pow(2.*M_PI,2.*3.-4.) );
}

void Decay2Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

namespace PHASIC {

  class Decay2_Channel_Generator: public Channel_Generator {
  public:
    
    Decay2_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new Decay2Channel(p_proc->NIn(),p_proc->NOut(),
				  (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }

  };// end of class Decay2_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(Decay2_Channel_Generator,"Decay2",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,Decay2_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new Decay2_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    Decay2_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"1->2 decay integrator";
}


NoChannel::NoChannel(int _nin,int _nout,Flavour * fl,Flavour res) :
  Single_Channel(_nin,_nout,fl)
{
  if (_nin != 2 || !(_nout==1 && fl[2].Kfcode()==kf_instanton)) {
    msg_Error()<<"Tried to initialize NoChannel for = "<<_nin<<" -> "<<_nout<<endl;
    Abort();
  }

  s      = smax  = pt2max = sqr(ATOOLS::rpa->gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  m_name   = "NoChannel";

  mass = width = 0.; 
  type = 0;
}

void NoChannel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran=0) {
  p[2] = p[0]+p[1];
}

void NoChannel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) { m_weight = 1.; }

void NoChannel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

std::string NoChannel::ChID() { return std::string("NoChannel"); }

namespace PHASIC {

  class No_Channel_Generator: public Channel_Generator {
  public:
    
    No_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new NoChannel(p_proc->NIn(),p_proc->NOut(),
			      (Flavour*)&p_proc->Flavours().front()));
      return 0;
    }

  };// end of class No_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(No_Channel_Generator,"NChannel",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,No_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new No_Channel_Generator(args);
}


void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    No_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"2->1 NoChannel integrator for Instanton production";
}
