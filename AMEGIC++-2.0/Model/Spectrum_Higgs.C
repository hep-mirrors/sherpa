#include "Spectrum_Higgs.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

void Spectrum_Higgs::Interface(Isajet* isa)
{
  CplEW.Init();
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->Higgs(_Z_R,_Z_H,v);

  msg.Tracking()<<"ISAJET v_Z:"<<v<<endl;  
    
  v = Flavour(kf::Z).Mass()*CplEW.SinThetaW()*CplEW.CosThetaW()/
    (2.*M_PI*M_PI*aqed->AqedFixed());

  msg.Tracking()<<"AMEGIC v_Z:"<<v<<endl;  

  tanb = _Z_H[0][0]/_Z_H[1][0];

  v1 = v*sqrt(1./(1.+sqr(tanb)));
  v2 = v1*tanb;  
  
  if (rpa.gen.Tracking())_Z_R.MatrixOut();
  msg.Tracking()<<"Z_H: "<<endl;  
  if (rpa.gen.Tracking())_Z_H.MatrixOut();
  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_Higgs::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());


  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  mu    = dr.GetValue<double>("mu");
  M_A0  = dr.GetValue<double>("M_A0");

  CplEW.Init();
  Masses_LO();
}

void Spectrum_Higgs::Masses_LO()
{
  double e       = sqrt(4.*M_PI*(aqed->AqedFixed()));
  double MW      = e/(2.*CplEW.SinThetaW())*v;
  double MZ      = e/(2.*CplEW.SinThetaW())*v/sqrt(1.-sqr(CplEW.SinThetaW()));

  v1 = v*sqrt(1./(1.+sqr(tanb)));
  v2 = v1*tanb;  

  // Charged Higgses
  double M_Hmin  = sqrt(sqr(MW)+sqr(M_A0));
  Flavour flav;
  flav = Flavour(kf::Hmin);flav.SetMass(M_Hmin);

  _Z_H[0][0] = v2/v;
  _Z_H[0][1] = -v1/v;
  _Z_H[1][0] = v1/v;
  _Z_H[1][1] = v2/v;
  //checked RK

  // Neutral Higgses
  Matrix<2> M;

  // hs = -mA^2*v1*v2/v^2


  M[0][0] = sqr(M_A0*v2/v)+sqr(MZ*v1/v);
  M[0][1] = -(sqr(M_A0)+sqr(MZ))*(v1*v2/sqr(v));
  M[1][0] = -(sqr(M_A0)+sqr(MZ))*(v1*v2/sqr(v));
  M[1][1] = sqr(M_A0*v1/v)+sqr(MZ*v2/v);
  //checked RK

  double evalues[2];

  M.DiagonalizeSort(evalues,_Z_R);

  flav = Flavour(kf::h0);flav.SetMass(sqrt(evalues[0]));
  flav = Flavour(kf::H0);flav.SetMass(sqrt(evalues[1]));
  flav = Flavour(kf::A0);flav.SetMass(M_A0);

  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"Higgs Masses :"<<endl;
  msg.Tracking()<<"m_h0 = "<<Flavour(kf::h0).Mass()<<endl;
  msg.Tracking()<<"m_H0 = "<<Flavour(kf::H0).Mass()<<endl;
  msg.Tracking()<<"m_A0 = "<<Flavour(kf::A0).Mass()<<endl;
  msg.Tracking()<<"m_H- = "<<Flavour(kf::Hmin).Mass()<<endl;
  msg.Tracking()<<"m_Z = "<<MZ<<endl;
  msg.Tracking()<<"v = "<<sqrt(sqr(v1)+sqr(v2))<<", tan(beta) = "<<v2/v1<<endl;
  msg.Tracking()<<"Test of tree--level relation : (m_h0^2+m_H0^2)/(m_A0^2+m_Z^2) = ";
  msg.Tracking()<<(evalues[0]+evalues[1])/(sqr(M_A0)+sqr(MZ))<<endl;
  msg.Tracking()<<"Z_R: "<<endl;  
  _Z_R.MatrixOut();
  msg.Tracking()<<"Z_H: "<<endl;  
  _Z_H.MatrixOut();
}





