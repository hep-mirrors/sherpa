#include "FSR_Channel.H"
#include "Channel_Elements.H"
#include "Channel_Basics.H"
#include "Run_Parameter.H"
#include "Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;


S1Channel::S1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize S1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = AMATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(AORGTOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "S-Channel 1";

  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void S1Channel::GeneratePoint(AMATOOLS::Vec4D * p,double * _ran=0) {
  CE.Isotropic2Momenta(p[0]+p[1],ms[2],ms[3],p[2],p[3],_ran[1],_ran[2]);
}

void S1Channel::GenerateWeight(AMATOOLS::Vec4D * p) {
  weight = 1. / ( CE.Isotropic2Weight(p[2],p[3]) * pow(2.*M_PI,2.*3.-4.) );
}

void S1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}





T1Channel::T1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize T1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = AMATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(AORGTOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "T-Channel 1";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void T1Channel::GeneratePoint(AMATOOLS::Vec4D * p,double * _ran =0) {
  CE.TChannelMomenta(p[0],p[1],p[2],p[3],ms[2],ms[3],0.,
		     0.5,0.,2.,1.,0,_ran[1],_ran[2]);
}

void T1Channel::GenerateWeight(AMATOOLS::Vec4D * p) {
  weight = 1. / ( CE.TChannelWeight(p[0],p[1],p[2],p[3],0.,0.5,0.,2.,1.,0) * pow(2.*M_PI,2*3.-4.) );
}

void T1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}





U1Channel::U1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize U1Channel with nout = "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = AMATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(AORGTOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "U-Channel 1";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void U1Channel::GeneratePoint(AMATOOLS::Vec4D * p,double * _ran =0) {
  CE.TChannelMomenta(p[0],p[1],p[3],p[2],ms[3],ms[2],0.,
		     0.5,0.,2.,1.,0,_ran[1],_ran[2]);
}

void U1Channel::GenerateWeight(AMATOOLS::Vec4D * p) {
  weight = 1. / ( CE.TChannelWeight(p[0],p[1],p[3],p[2],0.,0.5,0.,2.,1.,0) * pow(2.*M_PI,2*3.-4.) );
}

void U1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}




Decay2Channel::Decay2Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize Decay2Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = AMATOOLS::sqr(fl[i].Mass());
  rannum = 2;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(AORGTOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "Decay2-Channel 1";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void Decay2Channel::GeneratePoint(AMATOOLS::Vec4D * p,double * _ran=0) {
  CE.Isotropic2Momenta(p[0],ms[1],ms[2],p[1],p[2],_ran[0],_ran[1]);
}

void Decay2Channel::GenerateWeight(AMATOOLS::Vec4D * p) {
  weight = 1. / ( CE.Isotropic2Weight(p[1],p[2]) * pow(2.*M_PI,2.*3.-4.) );
}

void Decay2Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

