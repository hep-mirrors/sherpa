#include "Spectrum_sUpquarks.H"
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

void Spectrum_sUpquarks::Interface(Isajet* isa)
{
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->sUps(_Zu,_mu,_ws,_us);

  msg.Tracking()<<"Zu: "<<endl;
  if (rpa.gen.Tracking()) _Zu.MatrixOut();
  msg.Tracking()<<"ws: "<<endl;
  if (rpa.gen.Tracking()) _ws.MatrixOut();
  msg.Tracking()<<"us: "<<endl;
  if (rpa.gen.Tracking()) _us.MatrixOut();

  msg.Tracking()<<"======================================================"<<endl;

  isa->Gluino(mgluino);


  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_sUpquarks::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  mu    = dr.GetValue<double>("mu");
  mQ[0] = dr.GetValue<double>("MQ2(1)");
  mQ[1] = dr.GetValue<double>("MQ2(2)");
  mQ[2] = dr.GetValue<double>("MQ2(3)");
  mU[0] = dr.GetValue<double>("MU2(1)");
  mU[1] = dr.GetValue<double>("MU2(2)");
  mU[2] = dr.GetValue<double>("MU2(3)");
  m[0]  = dr.GetValue<double>("m_up");
  m[1]  = dr.GetValue<double>("m_charm");
  m[2]  = dr.GetValue<double>("m_top");
  A_t   = dr.GetValue<double>("A_t");

  Masses_LO();
}

void Spectrum_sUpquarks::Masses_LO()
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
  Matrix<3> D;
  Matrix<6> M;
  
  double helpA = -sqr(g/cos)*(sqr(v1)-sqr(v2))*(1.-4*sqr(cos))/24.;
  double helpB = sqr(e/cos)/6.*(sqr(v1)-sqr(v2));

  for (short int i=0;i<3;++i) {
    //double At = -(v1*ws[i]*ws[i]+v2*us[i]*us[i])/sqrt(2.);
    double MAt = 0.;
    if (i==2) MAt = A_t*m[i];
    for (short int j=0;j<3;++j) {
      A[i][j] = 0.;
      /*
      for (short int k=0;k<3;++k) {
	for (short int l=0;l<3;++l) {
	  A[i][j] += mQ[k]*mQ[l]*abs(CplEW.CKM(k,i)*conj(CplEW.CKM(l,j)));
	}
      }
      B[i][j] = mU[i]*mU[j];
      D[i][j] = -(v1*ws[i]*ws[j]+v2*us[i]*us[j])/sqrt(2.);
      */
      if (i==j) {
	for (short int k=0;k<3;++k) {
	  for (short int l=0;l<3;++l) {
	    A[i][i] += mQ[k]*mQ[l]*abs(CplEW.CKM(k,i)*conj(CplEW.CKM(l,j)));
	  }
	}
	B[i][i] = mU[i]*mU[i];
	D[i][i] = MAt;

	A[i][j] += helpA+sqr(m[i]);
	B[i][j] += helpB+sqr(m[i]);
	D[i][j] += -mu/tanb*m[i];
      }
    }
  }
  for (short int i=0;i<3;++i) {
    for (short int j=0;j<3;++j) {
      M[i][j]     = A[j][i];
      M[i+3][j+3] = B[i][j];
      M[i+3][j]   = D[i][j];
      M[i][j+3]   = D[j][i];
    }
  }

  double evalues[6];

  M.MatrixOut();
  M.Diagonalize(evalues,_Zu);
  Flavour flav;    
  flav = Flavour(kf::sUpL);flav.SetMass(sqrt(dabs(evalues[0])));
  flav = Flavour(kf::sUpR);flav.SetMass(sqrt(dabs(evalues[3])));
  flav = Flavour(kf::sCharmL);flav.SetMass(sqrt(dabs(evalues[1])));
  flav = Flavour(kf::sCharmR);flav.SetMass(sqrt(dabs(evalues[4])));
  flav = Flavour(kf::sTopL);flav.SetMass(sqrt(dabs(evalues[2])));
  flav = Flavour(kf::sTopR);flav.SetMass(sqrt(dabs(evalues[5])));

  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"sQuark masses :"<<endl;
  for (short int i=0;i<6;++i) msg.Tracking()<<"sUpquarks["<<i<<"] : "<<sqrt(dabs(evalues[i]))<<endl;
  _Zu.MatrixOut();
}








