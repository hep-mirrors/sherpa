#include "Spectrum_sDownquarks.H"
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


void Spectrum_sDownquarks::Interface(Isajet* isa)
{
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->sDowns(_Zd,_es,_ds);
  
  msg.Tracking()<<"Zd: "<<endl;
  if (rpa.gen.Tracking()) _Zd.MatrixOut();
  msg.Tracking()<<"es: "<<endl;
  if (rpa.gen.Tracking()) _es.MatrixOut();
  msg.Tracking()<<"ds: "<<endl;
  if (rpa.gen.Tracking()) _ds.MatrixOut();
  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_sDownquarks::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());


  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  mu    = dr.GetValue<double>("mu");
  mQ[0] = dr.GetValue<double>("MQ2(1)");
  mQ[1] = dr.GetValue<double>("MQ2(2)");
  mQ[2] = dr.GetValue<double>("MQ2(3)");
  mD[0] = dr.GetValue<double>("MD2(1)");
  mD[1] = dr.GetValue<double>("MD2(2)");
  mD[2] = dr.GetValue<double>("MD2(3)");
  m[0]  = dr.GetValue<double>("m_down");
  m[1]  = dr.GetValue<double>("m_strange");
  m[2]  = dr.GetValue<double>("m_bottom");
  A_b   = dr.GetValue<double>("A_b");

  Masses_LO();
}

void Spectrum_sDownquarks::Masses_LO()
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

  double helpA = -sqr(g/cos)*(sqr(v1)-sqr(v2))*(1.+2.*sqr(cos))/24.;
  double helpB = -sqr(e/cos)/12.*(sqr(v1)-sqr(v2));

  for (short int i=0;i<3;++i) {
    //double MAb = (v1*ds[i]*ds[j]-v2*es[i]*es[j])/sqrt(2.);
    double MAb = 0.;
    if (i==2) MAb = m[i]*A_b;
    for (short int j=0;j<3;++j) {
      /*
      A[i][j] = mQ[i]*mQ[j];
      B[i][j] = mD[i]*mD[j];
      C[i][j] = (v1*ds[i]*ds[j]-v2*es[i]*es[j])/sqrt(2.);
      */
      if (i==j) {
	A[i][j] = mQ[i]*mQ[j];
	B[i][j] = mD[i]*mD[j];
	C[i][j] = MAb;
	A[i][j] += helpA+sqr(m[i]);
	B[i][j] += helpB+sqr(m[i]);
	C[i][j] += -mu*tanb*m[i];
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

  M.Diagonalize(evalues,_Zd);

  //change evalues, that Left is always heavier than right
  
  if (evalues[2]<evalues[5]) {
    double help = evalues[2];
    evalues[2] = evalues[5];
    evalues[5] = help;
    for (short int i=0;i<6;i++) {
      help      = _Zd[i][2];
      _Zd[i][2] = _Zd[i][5];
      _Zd[i][5] = help;
    }
  }

  Flavour flav;    
  flav = Flavour(kf::sDownL);flav.SetMass(sqrt(dabs(evalues[0])));
  flav = Flavour(kf::sStrangeL);flav.SetMass(sqrt(dabs(evalues[1])));
  flav = Flavour(kf::sBottomL);flav.SetMass(sqrt(dabs(evalues[2])));
  flav = Flavour(kf::sDownR);flav.SetMass(sqrt(dabs(evalues[3])));
  flav = Flavour(kf::sStrangeR);flav.SetMass(sqrt(dabs(evalues[4])));
  flav = Flavour(kf::sBottomR);flav.SetMass(sqrt(dabs(evalues[5])));

  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"sQuark masses :"<<endl;
  for (short int i=0;i<6;++i) {
    msg.Tracking()<<"sDownquarks["<<i<<"] : ";
    msg.Tracking()<<sqrt(dabs(evalues[i]))<<endl;
  }
  _Zd.MatrixOut();
}








