#include "Run_Parameter.H"
#include "Message.H"
#include "Isajet.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace std;

extern "C" {
  void sugrun_(float&,float&,float&,float&,float&,float&);
  void chargino_(float&,float&,float&,float&,float&,float&);
  void neutralino_(float*,float*);
  void sneutrino_(float&,float&,float&);
  void higgses_(float&,float&,float&,float&,float&,float&,float&);
  void sups_(float*,float&,float&,float&);
  void sdowns_(float*,float&,float&);
  void sleptons_(float*,float&,float&,float&);
  void gluino_(float&);
};

Isajet::Isajet()
{
  //  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());
  Data_Read dr(rpa.GetPath()+std::string("/Isajet.dat"));

  if (!dr.GetValue<Switch::code>("Isajet")) return;
  float M0    = float(dr.GetValue<double>("M_0"));
  float M12   = float(dr.GetValue<double>("M_(1/2)"));
  float A0    = float(dr.GetValue<double>("A0"));
  float tanb  = float(dr.GetValue<double>("tanb"));
  float sgnmu = float(dr.GetValue<double>("sgn(mu)"));
  float Mt    = float(dr.GetValue<double>("M_t"));

  sugrun_(M0,M12,A0,tanb,sgnmu,Mt);  
}

void Isajet::Chargino(Matrix<2> &Zplus,Matrix<2> &Zminus)
{
  float MChi1,MChi2,gammaL,gammaR,ThX,ThY;
  
  chargino_(MChi1,MChi2,gammaL,gammaR,ThX,ThY);

  Flavour flav;
  flav = Flavour(kf::Chargino1);flav.SetMass(dabs(MChi1));
  if (MChi1<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Chargino2);flav.SetMass(dabs(MChi2));
  if (MChi2<0) flav.SetMassSign(-1);

  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"Charginomasses :"<<std::endl;
  msg.Tracking()<<"m_Chi_1 = "<<Flavour(kf::Chargino1).Mass()<<" Sign: "<<Flavour(kf::Chargino1).MassSign();
  msg.Tracking()<<", m_Chi_2 = "<<Flavour(kf::Chargino2).Mass()<<" Sign: "<<Flavour(kf::Chargino2).MassSign()<<std::endl;

  
  Zminus[0][0] = ::sin(gammaL);
  Zminus[1][0] = cos(gammaL);
  Zminus[0][1] = ThX*cos(gammaL);
  Zminus[1][1] = -ThX*::sin(gammaL);
  
  Zplus[0][0] = ::sin(gammaR);
  Zplus[1][0] = cos(gammaR);
  Zplus[0][1] = ThY*cos(gammaR);
  Zplus[1][1] = -ThY*::sin(gammaR);  
  

  //msg.Tracking()<<"ThX,THY: "<<ThX<<";"<<ThY<<std::endl;
  
}

void Isajet::Neutralino(Matrix<4> &ZN)
{
  float* M;
  M   = new float[4];
  float* Mix;
  Mix = new float[16];
  
  neutralino_(M,Mix);
  
  for (short int i=0;i<4;i++) M[i] = -M[i];
  //for (short int i=0;i<4;i++) M[i] = 100.;

  Matrix<4> M_Mix;

  for (short int i=0;i<4;i++) {
    for (short int j=0;j<4;j++) M_Mix[i][j] = Mix[i+j*4];
  }

  Flavour flav;
  flav = Flavour(kf::Neutralino1);flav.SetMass(dabs(M[0]));
  if (M[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Neutralino2);flav.SetMass(dabs(M[1]));
  if (M[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Neutralino3);flav.SetMass(dabs(M[2]));
  if (M[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Neutralino4);flav.SetMass(dabs(M[3]));
  if (M[3]<0) flav.SetMassSign(-1);

  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"Neutralinomasses :"<<std::endl;
  msg.Tracking()<<"m_Neu_1 = "<<Flavour(kf::Neutralino1).Mass();
  msg.Tracking()<<", m_Neu_2 = "<<Flavour(kf::Neutralino2).Mass()<<std::endl;
  msg.Tracking()<<"m_Neu_3 = "<<Flavour(kf::Neutralino3).Mass();
  msg.Tracking()<<", m_Neu_4 = "<<Flavour(kf::Neutralino4).Mass()<<std::endl;
  
  msg.Tracking()<<"Signs: "<<std::endl;
  msg.Tracking()<<"ms_Neu_1 = "<<Flavour(kf::Neutralino1).MassSign();
  msg.Tracking()<<", ms_Neu_2 = "<<Flavour(kf::Neutralino2).MassSign()<<std::endl;
  msg.Tracking()<<"ms_Neu_3 = "<<Flavour(kf::Neutralino3).MassSign();
  msg.Tracking()<<", ms_Neu_4 = "<<Flavour(kf::Neutralino4).MassSign()<<std::endl;

  for (short int i=0;i<4;i++) {
    //transposed with respect to Richardson !
    ZN[0][i] = M_Mix[3][i];
    ZN[1][i] = M_Mix[2][i];
    ZN[2][i] = -M_Mix[1][i];
    ZN[3][i] = -M_Mix[0][i];        
  }
  //for tests only
  //for (short int i=0;i<4;i++) {
  //  ZN[i][1] = - ZN[i][1];
  //  ZN[i][2] = - ZN[i][2];
  //}                                                                             

  //testing!!!!!!!
  //for (short int i=0;i<4;i++) {
  //  for (short int j=0;j<4;j++) ZN[i][j] = 1.; 
  //}

  delete[] M;
  delete[] Mix;
}

void Isajet::sNeutrino(Matrix<3> &ZNue)
{
  float mn1,mn2,mn3;
  sneutrino_(mn1,mn2,mn3);
  
  Flavour flav;    
  flav = Flavour(kf::sNu1);flav.SetMass(dabs(mn1));
  if (mn1<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sNu2);flav.SetMass(dabs(mn2));
  if (mn2<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sNu3);flav.SetMass(dabs(mn3));
  if (mn3<0) flav.SetMassSign(-1);
  
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ZNue[i][j] = 0.;
    ZNue[i][i] = 1.;
  }    

  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"sNeutrinomasses :"<<std::endl;
  msg.Tracking()<<"m_sNu_1 = "<<Flavour(kf::sNu1).Mass();
  msg.Tracking()<<", m_sNu_2 = "<<Flavour(kf::sNu2).Mass();
  msg.Tracking()<<", m_sNu_3 = "<<Flavour(kf::sNu3).Mass()<<std::endl;
}

void Isajet::Higgs(Matrix<2> &ZR,Matrix<2> & ZH,double& _v)
{
  float mh,mH,mA,mHC,tanb,alpha,vev;
  
  higgses_(mh,mH,mA,mHC,tanb,alpha,vev);

  _v = vev;
  
  double cosb = sqrt(1./(1.+sqr(tanb)));
  double sinb = cosb*tanb;  
  
  ZH[0][0] = sinb;
  ZH[0][1] = -cosb;
  ZH[1][0] = cosb;
  ZH[1][1] = sinb;
    
  double sina = ::sin(alpha);
  double cosa = cos(alpha);

  ZR[0][0] = -sina;
  ZR[0][1] = cosa;
  ZR[1][0] = cosa;
  ZR[1][1] = sina;
  
  Flavour flav;
  flav = Flavour(kf::h0);flav.SetMass(mh);
  flav = Flavour(kf::H0);flav.SetMass(mH);
  flav = Flavour(kf::A0);flav.SetMass(mA);
  flav = Flavour(kf::Hmin);flav.SetMass(mHC);

  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"Higgs Masses :"<<std::endl;
  msg.Tracking()<<"m_h0 = "<<Flavour(kf::h0).Mass()<<std::endl;
  msg.Tracking()<<"m_H0 = "<<Flavour(kf::H0).Mass()<<std::endl;
  msg.Tracking()<<"m_A0 = "<<Flavour(kf::A0).Mass()<<std::endl;
  msg.Tracking()<<"m_H- = "<<Flavour(kf::Hmin).Mass()<<std::endl;

  msg.Tracking()<<"Test of tree--level relation : (m_h0^2+m_H0^2)/(m_A0^2+m_Z^2) = ";
  msg.Tracking()<<(sqr(Flavour(kf::h0).Mass())+sqr(Flavour(kf::H0).Mass()))/
    (sqr(Flavour(kf::A0).Mass())+sqr(Flavour(kf::Z).Mass()))<<std::endl;
}

void Isajet::sUps(Matrix<6> &Zu,double &mu,Matrix<3> &ws,Matrix<3> &us)
{
  float* msups = new float[6];
  float thetaf,mu_1,_us22;
  sups_(msups,thetaf,mu_1,_us22);

  mu = mu_1;
  
  double us22  = _us22;
  
  double sinth = ::sin(thetaf);
  //msg.Tracking()<<"thUp :"<<thetaf<<std::endl;
  double costh = cos(thetaf);

  Flavour flav;    

  flav = Flavour(kf::sUpL);flav.SetMass(dabs(msups[0]));
  if (msups[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sCharmL);flav.SetMass(dabs(msups[1]));
  if (msups[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sTopL);flav.SetMass(dabs(msups[2]));
  if (msups[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sUpR);flav.SetMass(dabs(msups[3]));
  if (msups[3]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sCharmR);flav.SetMass(dabs(msups[4]));
  if (msups[4]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sTopR);flav.SetMass(dabs(msups[5]));
  if (msups[5]<0) flav.SetMassSign(-1);
  
  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"sQuark masses :"<<std::endl;
  for (short int i=0;i<6;++i) msg.Tracking()<<"sUpquarks["<<i<<"] : "<<dabs(msups[i])<<std::endl;


  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zu[i][j] = 0.;
    Zu[i][i] = 1.;
  }    
  Zu[2][2] = sinth;
  Zu[2][5] = costh;
  Zu[5][5] = sinth;
  Zu[5][2] = -costh;
  
  delete[] msups;

 for (short int i=0;i<2;i++) {
    for (short int j=0;j<2;j++) ws[i][j] = 0.;
   
  }    
  
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) us[i][j] = 0.;
    us[2][2] = us22;  
  }    
  
}

void Isajet::sDowns(Matrix<6> &Zd,Matrix<3> &es,Matrix<3> &ds)
{
  float* msdowns = new float[6];
  float thetaf,_ds22;
  sdowns_(msdowns,thetaf,_ds22);

  double ds22  = _ds22;
  double sinth = ::sin(thetaf);
  //msg.Tracking()<<"thDowns :"<<thetaf<<std::endl;
  double costh = cos(thetaf);

  Flavour flav;    

  flav = Flavour(kf::sDownL);flav.SetMass(dabs(msdowns[0]));
  if (msdowns[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sStrangeL);flav.SetMass(dabs(msdowns[1]));
  if (msdowns[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sBottomL);flav.SetMass(dabs(msdowns[2]));
  if (msdowns[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sDownR);flav.SetMass(dabs(msdowns[3]));
  if (msdowns[3]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sStrangeR);flav.SetMass(dabs(msdowns[4]));
  if (msdowns[4]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sBottomR);flav.SetMass(dabs(msdowns[5]));
  if (msdowns[5]<0) flav.SetMassSign(-1);
  
  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"sQuark masses :"<<std::endl;
  for (short int i=0;i<6;++i) msg.Tracking()<<"sDownquarks["<<i<<"] : "<<dabs(msdowns[i])<<std::endl;
  //msg.Tracking()<<"sQuark signs :"<<std::endl;
  //for (short int i=0;i<6;++i) msg.Tracking()<<"sDownquarks["<<i<<"] : "<<((msdowns[i] > 0) ? 1 : -1 )<<std::endl;

  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zd[i][j] = 0.;
    Zd[i][i] = 1.;
  }    

  Zd[2][2] = -sinth;
  Zd[2][5] = costh;
  Zd[5][5] = sinth;
  Zd[5][2] = costh;
  
  delete[] msdowns;

  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) es[i][j] = 0.;
  }    
  
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ds[i][j] = 0.;
    ds[2][2] = ds22;
  }    
}

void Isajet::sLeptons(Matrix<6> &Zl,double &mu,Matrix<3> &ls,Matrix<3> &ks)
{
  float* msleptons = new float[6];
  float thetaf,mu_1,_ls22;
  sleptons_(msleptons,thetaf,mu_1,_ls22);

  mu = mu_1;
 
  double ls22  = _ls22;
  double sinth = ::sin(thetaf);
  //msg.Tracking()<<"thLep :"<<thetaf<<std::endl;
  double costh = cos(thetaf);

  Flavour flav;    

  flav = Flavour(kf::sElectronL);flav.SetMass(dabs(msleptons[0]));
  if (msleptons[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sMuL);flav.SetMass(dabs(msleptons[1]));
  if (msleptons[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sTauL);flav.SetMass(dabs(msleptons[2]));
  if (msleptons[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sElectronR);flav.SetMass(dabs(msleptons[3]));
  if (msleptons[3]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sMuR);flav.SetMass(dabs(msleptons[4]));
  if (msleptons[4]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::sTauR);flav.SetMass(dabs(msleptons[5]));
  if (msleptons[5]<0) flav.SetMassSign(-1);
  
  msg.Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg.Tracking()<<"sLepton masses :"<<std::endl;
  for (short int i=0;i<6;++i) msg.Tracking()<<"sLeptons["<<i<<"] : "<<dabs(msleptons[i])<<std::endl;
   
  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zl[i][j] = 0.;
    Zl[i][i] = 1.;}    

  Zl[2][2] = sinth;
  Zl[2][5] = costh;
  Zl[5][5] = sinth;
  Zl[5][2] = -costh;
  
  delete[] msleptons;

  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ls[i][j] = 0.;}    
  ls[2][2] = ls22; 

  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ks[i][j] = 0.;
  }    
}

void Isajet::Gluino(double &m_GL)
{
  float _mGL;
  gluino_(_mGL);

  double mGL  = _mGL;
  
  Flavour flav;    

  flav = Flavour(kf::Gluino);flav.SetMass(dabs(mGL));
  if (mGL<0) flav.SetMassSign(-1);
  
  msg.Tracking()<<"Gluino mass : "<<dabs(mGL)<<std::endl;
        
}






