#include "FSR_Channel.H"
#include "Channel_Elements.H"
#include "Channel_Basics.H"
#include "Run_Parameter.H"
#include "Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


S1Channel::S1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize S1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "S-Channel";

  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(2,100,name,0);
}

void S1Channel::GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *cuts,double * _ran=0) {
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.Isotropic2Momenta(p[0]+p[1],ms[2],ms[3],p[2],p[3],ran[1],ran[2],-ctmax,ctmax);
}

void S1Channel::GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *cuts) {
  double ctmax=Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  double rans[2];
  weight = 1. / ( CE.Isotropic2Weight(p[2],p[3],rans[0],rans[1],-ctmax,ctmax) * pow(2.*M_PI,2.*3.-4.) );
  weight *= p_vegas->GenerateWeight(rans);
}

void S1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

std::string S1Channel::ChID() 
{
  return std::string("S-Channel");
}

T1Channel::T1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize T1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "T-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(2,100,name,0);
}

void T1Channel::GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *cuts,double * _ran =0) 
{
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.TChannelMomenta(p[0],p[1],p[2],p[3],ms[2],ms[3],0.,
		     .5,ctmax,-ctmax,1.,0,ran[1],ran[2]);
}

void T1Channel::GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *cuts) 
{
  double ctmax=Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  double rans[2];
  weight = 1. / ( CE.TChannelWeight(p[0],p[1],p[2],p[3],0.,
				    .5,ctmax,-ctmax,1.,0,rans[0],rans[1]) 
		  * pow(2.*M_PI,2*3.-4.) );
  weight *= p_vegas->GenerateWeight(rans);
}

void T1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}

std::string T1Channel::ChID() 
{
  return std::string("T-Channel");
}

U1Channel::U1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg.Error()<<"Tried to initialize U1Channel with nout = "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "U-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(2,100,name,0);
}

void U1Channel::GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *cuts,double * _ran =0) 
{
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=Min(cuts->cosmax[0][3],cuts->cosmax[1][2]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.TChannelMomenta(p[0],p[1],p[3],p[2],ms[3],ms[2],0.,
		     0.5,ctmax,-ctmax,1.,0,ran[1],ran[2]);
}

void U1Channel::GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data *cuts) 
{
  double ctmax=Min(cuts->cosmax[0][3],cuts->cosmax[1][2]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  double rans[2];
  weight = 1. / ( CE.TChannelWeight(p[0],p[1],p[3],p[2],0.,
				    .5,ctmax,-ctmax,1.,0,rans[0],rans[1]) 
		  * pow(2.*M_PI,2*3.-4.) );
  weight *= p_vegas->GenerateWeight(rans);
}

void U1Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}

std::string U1Channel::ChID() 
{
  return std::string("U-Channel");
}

Decay2Channel::Decay2Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=1) {
    msg.Error()<<"Tried to initialize Decay2Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "Decay2-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf::none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void Decay2Channel::GeneratePoint(ATOOLS::Vec4D * p,double * _ran=0) {
  CE.Isotropic2Momenta(p[0],ms[1],ms[2],p[1],p[2],_ran[1],_ran[2],-1.,1.);
}

void Decay2Channel::GenerateWeight(ATOOLS::Vec4D * p) {
  weight = 1. / ( CE.Isotropic2Weight(p[1],p[2],-1.,1.) * pow(2.*M_PI,2.*3.-4.) );
}

void Decay2Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

SimpleQCDChannel::SimpleQCDChannel(int _nin,int _nout,Flavour *flavs)
{  
  if (_nout!=2 || _nin!=2) abort();
  nin=_nin; 
  nout=_nout;
  ms = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i]=ATOOLS::sqr(flavs[i].Mass());
  rannum=3;
  rans = new double[rannum];
  name   = "S-Channel";
}

void SimpleQCDChannel::
GeneratePoint(ATOOLS::Vec4D *pi,ATOOLS::Cut_Data *cuts,double *ran) 
{
  Vec4D p=pi[0]+pi[1];
  double s1=ms[2], s2=ms[3];
  double s    = p.Abs2();
  double rs   = sqrt(dabs(s));
  Vec4D p1h;
  p1h[0]      = (s+s1-s2)/rs/2.;
  double p1m  = rs*Channel_Basics::SqLam(s,s1,s2)/2.;
  double phi  = 2.*M_PI*ran[2];
  double E1=(s+s1-s2)/2.0/sqrt(s);
  double sinthmin=cuts->etmin[2]/E1, sinthmax=sqrt(1.0-s1/sqr(E1));
  //  st=sinthmin*(1.-_ran[1])+_ran[1]*sinthmax;
  double st=CE.MasslessPropMomenta(1.0,sinthmin,sinthmax,ran[1]);
  double ct=sqrt(1.-st*st);
  Vec4D p1, p2;
  p1h = Vec4D(p1h[0],p1m*Vec3D(st*::sin(phi),st*cos(phi),ct));	
  Channel_Basics::Boost(0,p,p1h,p1);
  p2  = p+(-1.)*p1;
  pi[2]=p1;
  pi[3]=p2;
}

void SimpleQCDChannel::GenerateWeight(ATOOLS::Vec4D *p,ATOOLS::Cut_Data *cuts) 
{
  weight = 1. / ( CE.Isotropic2Weight(p[2],p[3],-1.,1.) * pow(2.*M_PI,2.*3.-4.) );
  double s1=ms[2], s2=ms[3];
  double s    = (p[0]+p[1]).Abs2();
  double E1=(s+s1-s2)/2.0/sqrt(s);
  double sinthmin=cuts->etmin[2]/E1, sinthmax=sqrt(1.0-s1/sqr(E1));
  //  weight/=p[2].CosTheta()/p[2].SinTheta()/(sinthmax-sinthmin);
  weight/=p[2].CosTheta()/p[2].SinTheta()*
    CE.MasslessPropWeight(1.0,sinthmin,sinthmax,p[2].SinTheta());
}

void SimpleQCDChannel::ISRInfo(int &type,double &mass,double &width) 
{
  type=2; 
  mass=sqrt(ms[2])+sqrt(ms[3]); 
  width=0.0;
}

std::string SimpleQCDChannel::ChID() 
{
  return "S-Channel";
}
