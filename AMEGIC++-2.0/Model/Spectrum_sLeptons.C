#include "Spectrum_sLeptons.H"
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

void Spectrum_sLeptons::Interface(Isajet* isa)
{
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->sLeptons(_ZLep,_mu,_ls,_ks);
  msg.Tracking()<<"ZLep: "<<endl;
  if (rpa.gen.Tracking()) _ZLep.MatrixOut();
  msg.Tracking()<<"ls: "<<endl;
  if (rpa.gen.Tracking()) _ls.MatrixOut();
  msg.Tracking()<<"ks: "<<endl;
  if (rpa.gen.Tracking()) _ks.MatrixOut();
  
  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_sLeptons::Init()
{
 Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  mu    = dr.GetValue<double>("mu");
  ml[0] = dr.GetValue<double>("ML2(1)");
  ml[1] = dr.GetValue<double>("ML2(2)");
  ml[2] = dr.GetValue<double>("ML2(3)");
  mr[0] = dr.GetValue<double>("MR2(1)");
  mr[1] = dr.GetValue<double>("MR2(2)");
  mr[2] = dr.GetValue<double>("MR2(3)");
  m[0]  = dr.GetValue<double>("m_e-");
  m[1]  = dr.GetValue<double>("m_mu");
  m[2]  = dr.GetValue<double>("m_tau");
  A_l   = dr.GetValue<double>("A_l");

  Masses_LO();
}

void Spectrum_sLeptons::Masses_LO()
{
  Couplings_EW CplEW;
  CplEW.Init();
  double e   = sqrt(4.*M_PI*(aqed->AqedFixed()));
  double g   = e/CplEW.SinThetaW();
  double cos = CplEW.CosThetaW();
  double sin = CplEW.SinThetaW();
  double v1 = sqrt(sqr(v)/(1.+sqr(tanb)));
  double v2 = v1*tanb;  

  Matrix<3> A;
  Matrix<3> B;
  Matrix<3> C;
  Matrix<6> M;

  double helpA = sqr(g/cos)*(sqr(v1)-sqr(v2))*(1.-2*sqr(cos))/8.;
  double helpB = -sqr(e/(2.*cos))*(sqr(v1)-sqr(v2));
  for (short int i=0;i<3;++i) {
    //double MAl = (v1*ls[i]*ls[j]-v2*ks[i]*ks[j])/sqrt(2.);
    double MAl = 0.;
    if (i==2) MAl = A_l*m[i];
    for (short int j=0;j<3;++j) {
      /*
      A[i][j] = ml[i]*ml[j];
      B[i][j] = mr[i]*mr[j];
      C[i][j] = (v1*ls[i]*ls[j]-v2*ks[i]*ks[j])/sqrt(2.);
      */
      if (i==j) {
	A[i][j] = ml[i]*ml[j];
	B[i][j] = mr[i]*mr[j];
	C[i][j] = MAl;
	A[i][j] += helpA+sqr(m[i]);
	B[i][j] += helpB+sqr(m[i]);
	C[i][j] += -tanb*mu*m[i];
      }
    }
  }
  for (short int i=0;i<3;++i) {
    for (short int j=0;j<3;++j) {
      M[i][j]     = A[j][i];
      M[i+3][j+3] = B[i][j];
      M[i+3][j]   = C[i][j];
      M[i][j+3]   = C[j][i];
    }
  }


  double evalues[6];

  //M.MatrixOut();
  M.Diagonalize(evalues,_ZLep);
  Flavour flav;    
  flav = Flavour(kf::sElectronL);flav.SetMass(sqrt(dabs(evalues[0])));
  flav = Flavour(kf::sElectronR);flav.SetMass(sqrt(dabs(evalues[3])));
  flav = Flavour(kf::sMuL);flav.SetMass(sqrt(dabs(evalues[1])));
  flav = Flavour(kf::sMuR);flav.SetMass(sqrt(dabs(evalues[4])));
  flav = Flavour(kf::sTauL);flav.SetMass(sqrt(dabs(evalues[2])));
  flav = Flavour(kf::sTauR);flav.SetMass(sqrt(dabs(evalues[5])));

  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"sLepton masses :"<<endl;
  for (short int i=0;i<6;++i) msg.Tracking()<<"SLepton["<<i<<"] : "<<sqrt(dabs(evalues[i]))<<endl;
  _ZLep.MatrixOut();
}


/*
void Spectrum_sLeptons::Isajet()
{
  string modelfilename=rpa.Get_Path()+string("/")+string("Isajet.dat");

  int parnum = 7;
  Data_Pointer** param = new Data_Pointer*[parnum];
  
  param[0] = new Data<Switch::code>;

  for (short int i=1;i<parnum;i++)
    param[i] = new Data<double>;

  param[0]->Set_Name(string("Isajet"));
  param[1]->Set_Name(string("m_selectron_L"));
  param[2]->Set_Name(string("m_smu_L"));
  param[3]->Set_Name(string("m_stau_L"));
  param[4]->Set_Name(string("m_selectron_R"));
  param[5]->Set_Name(string("m_smu_R"));
  param[6]->Set_Name(string("m_stau_R"));

  rpa.Read(modelfilename,param,parnum);

  Switch::code switch_tag;

  if (param[0]->Get_Value(switch_tag)==Switch::Off) return;

  double double_tag;

  Flavour flav;
  flav = Flavour(kf::sElectronL);flav.SetMass(param[1]->Get_Value(double_tag));
  flav = Flavour(kf::sElectronR);flav.SetMass(param[4]->Get_Value(double_tag));
  flav = Flavour(kf::sMuL);flav.SetMass(param[2]->Get_Value(double_tag));
  flav = Flavour(kf::sMuR);flav.SetMass(param[5]->Get_Value(double_tag));
  flav = Flavour(kf::sTauL);flav.SetMass(param[3]->Get_Value(double_tag));
  flav = Flavour(kf::sTauR);flav.SetMass(param[6]->Get_Value(double_tag));


  for (short int i=0;i<6;++i) {
    for (short int j=0;j<6;++j) _ZLep[i][j] = 0.;
    _ZLep[i][i] = 1.;
  }

  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"sLepton masses :"<<endl;
  for (short int i=1;i<7;++i) msg.Tracking()<<"SLepton["<<i<<"] : "<<param[i]->Get_Value(double_tag)<<endl;
  msg.Tracking()<<"Mass matrix diagonal ... no slepton mixing !"<<endl;

  _ZLep.MatrixOut();
}
*/






