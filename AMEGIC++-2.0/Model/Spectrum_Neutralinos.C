#include "Spectrum_Neutralinos.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

void Spectrum_Neutralinos::Interface(Isajet* isa)
{
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->Neutralino(_Z_N);

  msg.Tracking()<<"ZNeutralino: "<<endl;
  if (rpa.gen.Tracking()) _Z_N.MatrixOut();
  msg.Tracking()<<"ZN(trans): "<<endl;
  if (rpa.gen.Tracking()) (_Z_N.Dagger()).MatrixOut();
  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_Neutralinos::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());


  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  m_2   = dr.GetValue<double>("m_2");
  m_3   = dr.GetValue<double>("m_3");
  mu    = dr.GetValue<double>("mu");

  CplEW.Init();
  Masses_LO();
}

void Spectrum_Neutralinos::Masses_LO()
{
  double e   = sqrt(4.*M_PI*(aqed->AqedFixed()));

  double v1  = v*sqrt(1./(1.+sqr(tanb)));
  double v2  = v1*tanb;  
  double sin = CplEW.SinThetaW();
  double cos = sqrt(1.-sqr(sin));
  Matrix<4> M;

  M[0][0] = 2.*m_3;
  M[0][1] = 0.;
  M[0][2] = -e*v1/(2.*cos);
  M[0][3] = e*v2/(2.*cos);
  M[1][0] = 0.;
  M[1][1] = 2.*m_2;
  M[1][2] = e*v1/(2.*sin);
  M[1][3] = -e*v2/(2.*sin);
  M[2][0] = -e*v1/(2.*cos);
  M[2][1] = e*v1/(2.*sin);
  M[2][2] = 0.;
  M[2][3] = -mu;
  M[3][0] = e*v2/(2.*cos);
  M[3][1] = -e*v2/(2.*sin);
  M[3][2] = -mu;
  M[3][3] = 0.;


  double evalues[4];

  M.DiagonalizeSort(evalues,_Z_N);

  Flavour flav;
  flav = Flavour(kf::Neutralino1);flav.SetMass(dabs(evalues[0]));
  if (evalues[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Neutralino2);flav.SetMass(dabs(evalues[1]));
  if (evalues[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Neutralino3);flav.SetMass(dabs(evalues[2]));
  if (evalues[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Neutralino4);flav.SetMass(dabs(evalues[3]));
  if (evalues[3]<0) flav.SetMassSign(-1);

  //Matrix<4> help = _Z_N.Dagger()*M*_Z_N;
  _Z_N.MatrixOut();
  msg.Tracking()<<"Signs: "<<endl;
  msg.Tracking()<<"ms_Neu_1 = "<<Flavour(kf::Neutralino1).MassSign();
  msg.Tracking()<<", ms_Neu_2 = "<<Flavour(kf::Neutralino2).MassSign()<<endl;
  msg.Tracking()<<"ms_Neu_3 = "<<Flavour(kf::Neutralino3).MassSign();
  msg.Tracking()<<", ms_Neu_4 = "<<Flavour(kf::Neutralino4).MassSign()<<endl;
  


  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"Neutralinomasses :"<<endl;
  msg.Tracking()<<"m_Neu_1 = "<<Flavour(kf::Neutralino1).Mass();
  msg.Tracking()<<", m_Neu_2 = "<<Flavour(kf::Neutralino2).Mass()<<endl;
  msg.Tracking()<<"m_Neu_3 = "<<Flavour(kf::Neutralino3).Mass();
  msg.Tracking()<<", m_Neu_4 = "<<Flavour(kf::Neutralino4).Mass()<<endl;
}




