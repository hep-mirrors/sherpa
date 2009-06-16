#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Channel_Basics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

#include <stdio.h>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


S1Channel::S1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg_Error()<<"Tried to initialize S1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 2;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "S-Channel";

  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(rannum,100,name,0);
}

S1Channel::~S1Channel()
{
  delete p_vegas;
}

void S1Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran=0) {
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.Isotropic2Momenta(p[0]+p[1],ms[2],ms[3],p[2],p[3],ran[0],ran[1],-ctmax,ctmax);
}

void S1Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) {
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
    msg_Error()<<"Tried to initialize T1Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; 
  nout = _nout;
  ms   = new double[nin+nout];
  for (int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3*nout-4;
  rans   = new double[rannum];
  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.0;
  E      = 0.5 * sqrt(s);
  name   = "T-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(rannum,100,name,0);
}

T1Channel::~T1Channel()
{
  delete p_vegas;
}

void T1Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran =0) 
{
  double ctmax=Min(cuts->cosmax[0][2],cuts->cosmax[1][3]);
  double *ran = p_vegas->GeneratePoint(_ran);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.TChannelMomenta(p[0],p[1],p[2],p[3],ms[2],ms[3],0.,
		     .5,ctmax,-ctmax,1.,0,ran[0],ran[1]);
}

void T1Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) 
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
  return name;
}

T2Channel::T2Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  nin  = _nin; 
  nout = _nout;
  ms   = new double[nin+nout];
  for (int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3*nout-4;
  rans   = new double[rannum];
  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.0;
  E      = 0.5 * sqrt(s);
  name   = "T2-Channel";
  if (nout>2) name   = ToString(nin)+"->"+ToString(nout)+"_"+name;
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(rannum,100,name,0);
}

T2Channel::~T2Channel()
{
  delete p_vegas;
}

void T2Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran =0) 
{
  double ctmax=1.0;
  double *ran = p_vegas->GeneratePoint(_ran);
  Vec4D pa, pr(p[1]);
  double tmass=0.0, rsmax=sqrt((p[0]+p[1]).Abs2()) ;
  std::string tagl;
  for (int i=nin;i<nin+nout-1;++i) tagl=tagl+ToString(i);
  for (int i=0;i<nout-2;++i) {
    double ctmin=-(ctmax=1.0);
    double smin = cuts->Getscut(tagl);
    double smax = sqr(rsmax-sqrt(ms[nin+nout-i-1]));
    double s = CE.MasslessPropMomenta(.5,smin,smax,ran[3*i]);
    CE.TChannelMomenta(p[0],pr,pa,p[nin+nout-i-1],s,ms[nin+nout-i-1],tmass,
		       .5,ctmax,ctmin,1.,0,ran[3*i+1],ran[3*i+2]);
    pr-=p[nin+nout-i-1];
    rsmax=sqrt((p[0]+pr).Abs2());
    tagl=tagl.substr(0,tagl.length()-1);
  }
  CE.TChannelMomenta(p[0],pr,p[2],p[3],ms[2],ms[3],tmass,
		     .5,ctmax,-ctmax,1.,0,ran[3*nout-4-2],ran[3*nout-4-1]);
}

void T2Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) 
{
  double ctmax=1.0;
  double wt = 1.0;
  Vec4D pa(p[0]+p[1]), pr(p[1]);
  double tmass=0.0, rsmax=sqrt((p[0]+p[1]).Abs2()) ;
  std::string tagl;
  for (int i=nin;i<nin+nout-1;++i) tagl=tagl+ToString(i);
  for (int i=0;i<nout-2;++i) {
    double ctmin=-(ctmax=1.0);
    double smin = cuts->Getscut(tagl);
    double smax = sqr(rsmax-sqrt(ms[nin+nout-i-1]));
    pa-=p[nin+nout-i-1];
    wt *= CE.MasslessPropWeight(.5,smin,smax,pa.Abs2(),rans[3*i]);
    wt *= CE.TChannelWeight(p[0],pr,pa,p[nin+nout-i-1],tmass,
			    .5,ctmax,ctmin,1.,0,rans[3*i+1],rans[3*i+2]);
    pr-=p[nin+nout-i-1];
    rsmax=sqrt((p[0]+pr).Abs2());
    tagl=tagl.substr(0,tagl.length()-1);
  }
  wt *= CE.TChannelWeight(p[0],pr,p[2],p[3],tmass,.5,ctmax,-ctmax,1.,0,
			  rans[3*nout-4-2],rans[3*nout-4-1]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,nout*3.-4.);
  weight = wt;
  weight *= p_vegas->GenerateWeight(rans);
}

void T2Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}

std::string T2Channel::ChID() 
{
  return name;
}

T3Channel::T3Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  nin  = _nin; 
  nout = _nout;
  ms   = new double[nin+nout];
  for (int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 3*nout-4;
  rans   = new double[rannum];
  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.0;
  E      = 0.5 * sqrt(s);
  name   = "T3-Channel";
  if (nout>2) name   = ToString(nin)+"->"+ToString(nout)+"_"+name;
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(rannum,100,name,0);
}

T3Channel::~T3Channel()
{
  delete p_vegas;
}

void T3Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran =0) 
{
  double ctmax=1.0;
  double *ran = p_vegas->GeneratePoint(_ran);
  Vec4D pa, pl(p[0]);
  double tmass=0.0, rsmax=sqrt((p[0]+p[1]).Abs2()) ;
  std::string tagl;
  for (int i=nin;i<nin+nout-1;++i) tagl=tagl+ToString(i);
  for (int i=0;i<nout-2;++i) {
    double ctmin=-(ctmax=1.0);
    double smin = cuts->Getscut(tagl);
    double smax = sqr(rsmax-sqrt(ms[nin+i]));
    double s = CE.MasslessPropMomenta(.5,smin,smax,ran[3*i]);
    CE.TChannelMomenta(pl,p[1],p[nin+i],pa,ms[nin+i],s,tmass,
		       .5,ctmax,ctmin,1.,0,ran[3*i+1],ran[3*i+2]);
    pl-=p[nin+i];
    rsmax=sqrt((pl+p[1]).Abs2());
    tagl=tagl.substr(1);
  }
  CE.TChannelMomenta(pl,p[1],p[nin+nout-2],p[nin+nout-1],
		     ms[nin+nout-2],ms[nin+nout-1],tmass,
		     .5,ctmax,-ctmax,1.,0,ran[3*nout-4-2],ran[3*nout-4-1]);
}

void T3Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) 
{
  double ctmax=1.0;
  double wt = 1.0;
  Vec4D pa(p[0]+p[1]), pl(p[0]);
  double tmass=0.0, rsmax=sqrt((p[0]+p[1]).Abs2()) ;
  std::string tagl;
  for (int i=nin;i<nin+nout-1;++i) tagl=tagl+ToString(i);
  for (int i=0;i<nout-2;++i) {
    double ctmin=-(ctmax=1.0);
    double smin = cuts->Getscut(tagl);
    double smax = sqr(rsmax-sqrt(ms[nin+nout-i-1]));
    pa-=p[nin+i];
    wt *= CE.MasslessPropWeight(.5,smin,smax,pa.Abs2(),rans[3*i]);
    wt *= CE.TChannelWeight(pl,p[1],pa,p[nin+i],tmass,
			    .5,ctmax,ctmin,1.,0,rans[3*i+1],rans[3*i+2]);
    pl-=p[nin+i];
    rsmax=sqrt((pl+p[1]).Abs2());
    tagl=tagl.substr(1);
  }
  wt *= CE.TChannelWeight(pl,p[1],p[nin+nout-2],p[nin+nout-1],
			  tmass,.5,ctmax,-ctmax,1.,0,
			  rans[3*nout-4-2],rans[3*nout-4-1]);
  if (wt!=0.) wt = 1.0/wt/pow(2.*M_PI,nout*3.-4.);
  weight = wt;
  weight *= p_vegas->GenerateWeight(rans);
}

void T3Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = 0; _mass = mass; _width = width;
}

std::string T3Channel::ChID() 
{
  return name;
}

U1Channel::U1Channel(int _nin,int _nout,Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=2) {
    msg_Error()<<"Tried to initialize U1Channel with nout = "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 2;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "U-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
  p_vegas = new Vegas(rannum,100,name,0);
}

U1Channel::~U1Channel()
{
  delete p_vegas;
}

void U1Channel::GeneratePoint(ATOOLS::Vec4D * p,Cut_Data *cuts,double * _ran =0) 
{
  double *ran = p_vegas->GeneratePoint(_ran);
  double ctmax=Min(cuts->cosmax[0][3],cuts->cosmax[1][2]);
  double s=(p[0]+p[1]).Abs2(), E12=sqr(s+ms[2]-ms[3])/4.0/s;
  ctmax=Min(ctmax,sqrt(1.0-sqr(cuts->etmin[2])/E12));
  CE.TChannelMomenta(p[0],p[1],p[3],p[2],ms[3],ms[2],0.,
		     0.5,ctmax,-ctmax,1.,0,ran[0],ran[1]);
}

void U1Channel::GenerateWeight(ATOOLS::Vec4D * p,Cut_Data *cuts) 
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

Decay2Channel::Decay2Channel(int _nin,int _nout,const Flavour * fl,Flavour res) 
{  
  if (_nout != 2 || _nin!=1) {
    msg_Error()<<"Tried to initialize Decay2Channel with nout = "<<_nin<<" -> "<<_nout<<endl;
    abort();
  }
  nin  = _nin; nout = _nout;
  ms   = new double[nin+nout];
  for (short int i=0;i<nin+nout;i++) ms[i] = ATOOLS::sqr(fl[i].Mass());
  rannum = 2;
  rans   = new double[rannum];

  s      = smax  = pt2max = sqr(ATOOLS::rpa.gen.Ecms());
  pt2min = 0.;
  E      = 0.5 * sqrt(s);
  name   = "Decay2-Channel";
  mass   = width = 0.; 
  type   = 0;
  if (res!=Flavour(kf_none)) {
    mass = res.Mass(); width = res.Width(); type = 1;
  }
}

void Decay2Channel::GeneratePoint(ATOOLS::Vec4D * p,double * _ran=0) {
  CE.Isotropic2Momenta(p[0],ms[1],ms[2],p[1],p[2],_ran[0],_ran[1],-1.,1.);
}

void Decay2Channel::GenerateWeight(ATOOLS::Vec4D * p) {
  weight = 1. / ( CE.Isotropic2Weight(p[1],p[2],-1.,1.) * pow(2.*M_PI,2.*3.-4.) );
}

void Decay2Channel::ISRInfo(int & _type,double & _mass,double & _width) {
  _type = type; _mass = mass; _width = width;
}

