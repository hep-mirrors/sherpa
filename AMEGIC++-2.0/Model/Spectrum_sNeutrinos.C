#include "Spectrum_sNeutrinos.H"
#include "Couplings_EW.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

void Spectrum_sNeutrinos::Interface(Isajet* isa)
{
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->sNeutrino(_Znue);
  
  msg.Tracking()<<"ZsNeutrino: "<<endl;
  if (rpa.gen.Tracking()) _Znue.matrix_out();
  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_sNeutrinos::Init()
{
 Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  ml[0] = dr.GetValue<double>("ML2(1)");
  ml[1] = dr.GetValue<double>("ML2(2)");
  ml[2] = dr.GetValue<double>("ML2(3)");

  Masses_LO();
}

void Spectrum_sNeutrinos::Masses_LO()
{
  Couplings_EW CplEW;
  CplEW.Init();
  double e  = sqrt(4.*M_PI*(aqed->AqedFixed()));

  double v1 = sqrt(sqr(v)/(1.+sqr(tanb)));
  double v2 = v1*tanb;  

  Matrix<3> M;

  //mL2 Matrix = Zero
  Matrix<3> mL2;

  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) mL2[i][j] = ml[i]*ml[j];
  }

  //diagonal -- test
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) mL2[i][j] = 0.;
    mL2[i][i] = ml[i]*ml[i];
  }


  //msg.Tracking()<<"ML - Matrix: "<<endl;
  //mL2.matrix_out();
 

  for (short int i=0;i<3;i++) {
    M[i][i] = sqr(e)*(sqr(v1)-sqr(v2))/
      (8.*sqr(CplEW.SinThetaW())*sqr(CplEW.CosThetaW()));
    for (short int j=0;j<3;j++) M[i][j] += mL2[i][j];
  }
  double evalues[3];

  M.Diagonalize(evalues,_Znue);

  Flavour flav;    
  flav = Flavour(kf::sNu1);flav.set_mass(sqrt(dabs(evalues[0])));
  flav = Flavour(kf::sNu2);flav.set_mass(sqrt(dabs(evalues[1])));
  flav = Flavour(kf::sNu3);flav.set_mass(sqrt(dabs(evalues[2])));

  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"sNeutrinomasses :"<<endl;
  msg.Tracking()<<"m_sNu_1 = "<<Flavour(kf::sNu1).mass();
  msg.Tracking()<<", m_sNu_2 = "<<Flavour(kf::sNu2).mass();
  msg.Tracking()<<", m_sNu_3 = "<<Flavour(kf::sNu3).mass()<<endl;
  _Znue.matrix_out();
}
