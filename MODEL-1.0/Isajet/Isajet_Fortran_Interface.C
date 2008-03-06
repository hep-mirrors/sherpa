#include "Isajet_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MathTools.H"

using namespace ISAJET;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

extern "C" {
  void isasusyinter_(float *,float *,float *,float *,int &,const char *);
  void chargino_(float&,float&,float&,float&,float&,float&);
  void neutralino_(float*,float*);
  void sneutrino_(float&,float&,float&);
  void higgses_(float&,float&,float&,float&,float&,float&,float&);
  void sups_(float*,float&,float&,float&);
  void sdowns_(float*,float&,float&);
  void sleptons_(float*,float&,float&,float&);
  void gluino_(float&);
  void filldecays_(int &,float &,int &,int *,int *,float *);
}




Isajet_Fortran_Interface::Isajet_Fortran_Interface(Data_Read * _dataread,
						   Model_Base * _model) :
  Spectrum_Generator_Base(_dataread,_model) { }



Isajet_Fortran_Interface::~Isajet_Fortran_Interface() { }



void Isajet_Fortran_Interface::Run(std::string _model) {
  float * xxsm = new float[15];
  xxsm[0]   = float(p_model->ScalarConstant(std::string("Yukawa_u")));
  xxsm[1]   = float(p_model->ScalarConstant(std::string("Yukawa_d")));
  xxsm[2]   = float(p_model->ScalarConstant(std::string("Yukawa_s")));
  xxsm[3]   = float(p_model->ScalarConstant(std::string("Yukawa_c")));
  xxsm[4]   = float(p_model->ScalarConstant(std::string("Yukawa_b")));
  xxsm[5]   = float(p_model->ScalarConstant(std::string("Yukawa_t")));
  xxsm[6]   = float(p_model->ScalarConstant(std::string("Yukawa_e")));
  xxsm[7]   = float(p_model->ScalarConstant(std::string("Yukawa_mu")));
  xxsm[8]   = float(p_model->ScalarConstant(std::string("Yukawa_tau")));
  xxsm[9]   = float(p_model->ScalarConstant(std::string("MZ")));
  xxsm[10]  = float(p_model->ScalarConstant(std::string("GammaW")));
  xxsm[11]  = float(p_model->ScalarConstant(std::string("GammaZ")));
  xxsm[12]  = float(p_model->ScalarFunction(std::string("alpha_QED"),sqr(xxsm[9]))); 
  xxsm[13]  = float(p_model->ScalarConstant(std::string("sin2_thetaW"))); 
  xxsm[14]  = float(p_model->ScalarFunction(std::string("alpha_S"),sqr(xxsm[9]))); 

  float * xsugin = new float[7];
  for (int i=0;i<7;i++)    xsugin[i]=0.;
  float * xgmin = new float[14];
  for (int i=0;i<14;i++)   xgmin[i]=0.;
  float * xnusug = new float[18];
  for (int i=0;i<18;i++)   xnusug[i]=1.e20;
  
  int model;
  if (_model==std::string("mSUGRA") ||
      _model==std::string("non-universal SUGRA") ||
      _model==std::string("SUGRA with enforced unification") ) {
    model     = 1;
    xsugin[0] = float(p_model->ScalarConstant(std::string("m0")));
    xsugin[1] = float(p_model->ScalarConstant(std::string("m12")));
    xsugin[2] = float(p_model->ScalarConstant(std::string("A0")));
    xsugin[3] = float(p_model->ScalarConstant(std::string("tan(beta)")));
    xsugin[4] = float(p_model->ScalarNumber(std::string("sign(mu)")));
    xsugin[5] = float(p_model->ScalarConstant(std::string("mT")));
    if (_model==std::string("non-universal SUGRA")) {
      model     = 3;
      xsugin[6] = float(p_model->ScalarConstant(std::string("m_SUSY")));

      xnusug[0]  = float(p_model->ScalarConstant(std::string("m_U1")));
      xnusug[1]  = float(p_model->ScalarConstant(std::string("m_SU2")));
      xnusug[2]  = float(p_model->ScalarConstant(std::string("m_SU3")));
      xnusug[3]  = float(p_model->ScalarConstant(std::string("A(tau)_GUT")));
      xnusug[4]  = float(p_model->ScalarConstant(std::string("A(b)_GUT")));
      xnusug[5]  = float(p_model->ScalarConstant(std::string("A(t)_GUT")));
      xnusug[6]  = float(p_model->ScalarConstant(std::string("H(d)_GUT")));
      xnusug[7]  = float(p_model->ScalarConstant(std::string("H(u)_GUT")));
      xnusug[8]  = float(p_model->ScalarConstant(std::string("M(e1R)_GUT")));
      xnusug[9]  = float(p_model->ScalarConstant(std::string("M(l1)_GUT")));
      xnusug[10] = float(p_model->ScalarConstant(std::string("M(d1R)_GUT")));
      xnusug[11] = float(p_model->ScalarConstant(std::string("M(u1R)_GUT")));
      xnusug[12] = float(p_model->ScalarConstant(std::string("M(q1)_GUT")));
      xnusug[13]  = float(p_model->ScalarConstant(std::string("M(e3R)_GUT")));
      xnusug[14]  = float(p_model->ScalarConstant(std::string("M(l3)_GUT")));
      xnusug[15] = float(p_model->ScalarConstant(std::string("M(d3R)_GUT")));
      xnusug[16] = float(p_model->ScalarConstant(std::string("M(u3R)_GUT")));
      xnusug[17] = float(p_model->ScalarConstant(std::string("M(q3)_GUT")));
    }
    if(_model==std::string("SUGRA with enforced unification")) {
      model     = 4;
    }
  }
  if (_model==std::string("AMSB")) {
    model     = 7;
    xsugin[0] = float(p_model->ScalarConstant(std::string("m0")));
    xsugin[1] = float(p_model->ScalarConstant(std::string("m32")));
    xsugin[2] = 0.;
    xsugin[3] = float(p_model->ScalarConstant(std::string("tan(beta)")));
    xsugin[4] = float(p_model->ScalarNumber(std::string("sign(mu)")));
    xsugin[5] = float(p_model->ScalarConstant(std::string("mT")));
  }
  if (_model==std::string("mGMSB") ||
      _model==std::string("non-minimal GMSB")) {
    model      = 2;
    xsugin[0]  = float(p_model->ScalarConstant(std::string("Lambda_m")));
    xsugin[1]  = float(p_model->ScalarConstant(std::string("m_mes")));
    xsugin[2]  = float(p_model->ScalarNumber(std::string("n_mes")));
    xsugin[3]  = float(p_model->ScalarConstant(std::string("tan(beta)")));
    xsugin[4]  = float(p_model->ScalarNumber(std::string("sign(mu)")));
    xsugin[5]  = float(p_model->ScalarConstant(std::string("mT")));
    xsugin[6]  = float(p_model->ScalarConstant(std::string("c_grav")));
    if(_model==std::string("non-minimal GMSB")) {
      model     = 5;
      xgmin[7]  = float(p_model->ScalarConstant(std::string("c_gauge")));
      xgmin[8]  = float(p_model->ScalarConstant(std::string("Delta_Hd")));
      xgmin[9]  = float(p_model->ScalarConstant(std::string("Delta_Hu")));
      xgmin[10] = float(p_model->ScalarConstant(std::string("Delta_Y")));
      xgmin[11] = float(p_model->ScalarConstant(std::string("n5_1")));
      xgmin[12] = float(p_model->ScalarConstant(std::string("n5_2")));
      xgmin[13] = float(p_model->ScalarConstant(std::string("n5_3")));
    }
  }
  
  const char * help;
  std::string full = p_dataread->GetValue<std::string>("OUTPUTFILE",std::string("Isajet.out"));
  help = full.c_str();
  isasusyinter_(xsugin,xnusug,xgmin,xxsm,model,help);

  delete [] xxsm;
  delete [] xgmin;
  delete [] xsugin;
  delete [] xnusug;
}


void Isajet_Fortran_Interface::FillMasses() {
  CharginoMasses();
  NeutralinoMasses();
  GluinoMasses();
  HiggsMasses();
  sNeutrinoMasses();
  sUpMasses();
  sDownMasses();
  sLeptonMasses();
}

void Isajet_Fortran_Interface::FillDecays() {
  msg_Info()<<"---------- Set particles widths according to Isasusy ! ---------"<<endl;
  
  Decays(Flavour(kf_Chargino1));
  Decays(Flavour(kf_Chargino2));
  
  Decays(Flavour(kf_Neutralino1));
  Decays(Flavour(kf_Neutralino2));
  Decays(Flavour(kf_Neutralino3));
  Decays(Flavour(kf_Neutralino4));

  Decays(Flavour(kf_Gluino));
  
  Decays(Flavour(kf_h0));
  Decays(Flavour(kf_H0));
  Decays(Flavour(kf_A0));
  Decays(Flavour(kf_Hmin));

  Decays(Flavour(kf_sElectronL));
  Decays(Flavour(kf_sMuL));
  Decays(Flavour(kf_sTau2));
  Decays(Flavour(kf_sElectronR));
  Decays(Flavour(kf_sMuR));
  Decays(Flavour(kf_sTau1));

  Decays(Flavour(kf_sNu1));
  Decays(Flavour(kf_sNu2));
  Decays(Flavour(kf_sNu3));
  
  Decays(Flavour(kf_sUpL));
  Decays(Flavour(kf_sCharmL));
  Decays(Flavour(kf_sTop2));
  Decays(Flavour(kf_sUpR));
  Decays(Flavour(kf_sCharmR));
  Decays(Flavour(kf_sTop1));
  
  Decays(Flavour(kf_sDownL));
  Decays(Flavour(kf_sStrangeL));
  Decays(Flavour(kf_sBottom2));
  Decays(Flavour(kf_sDownR));
  Decays(Flavour(kf_sStrangeR));
  Decays(Flavour(kf_sBottom1));
}

int Isajet_Fortran_Interface::FlavourToIsaID(Flavour flav) {

  int isaID = 0;

  int sherpaID = int(flav);
  bool antifl = 0;

  if (sherpaID<0) {
    sherpaID*=-1;
    antifl   = 1;
  }
  
  switch(sherpaID) {
  case 41 : isaID = 39; //39 equals chargino^-_1
    break;
  case 42 : isaID = 49; //49 equals chargino^-_2
    break;
  case 43 : isaID = 30; //30 equals neutralino_1
    break;
  case 44 : isaID = 40; //40 equals neutralino_2
    break;
  case 45 : isaID = 50; //50 equals neutralino_3
    break;
  case 46 : isaID = 60; //60 equals neutralino_4
    break;
  case 47 : isaID = 29; //29 equals gluino
    break;
  case 31 : isaID = 82; //82 equals h0
    break;  
  case 32 : isaID = 83; //83 equals H0
    break;
  case 33 : isaID = 84; //84 equals A0
    break;
  case 34 : isaID = 86; //86 equals anti-H^-
    break;
  case 81 : isaID = 31; //31 equals sNu1
     break;
  case 82 : isaID = 33; //33 equals sNu2
    break;
  case 83 : isaID = 35; //35 equals sNu3
    break;
  case 71 : isaID = 32; //32 equals sElectronL
    break;
  case 72 : isaID = 34; //34 equals sMuonL
    break;  
  case 73 : isaID = 56; //36 equals sTau2
    break;  
  case 74 : isaID = 52; //52 equals sElectronR
    break;  
  case 75 : isaID = 54; //54 equals sMuonR
    break;  
  case 76 : isaID = 36; //56 equals sTau1
    break;
  case 51 : isaID = 21; //21 equals sUpL
    break;
  case 52 : isaID = 24; //24 equals sCharmL
    break;
  case 53 : isaID = 46; //26 equals sTop2
    break;
  case 54 : isaID = 41; //41 equals sUpR
    break;
  case 55 : isaID = 44; //44 equals sCharmR
    break;
  case 56 : isaID = 26; //46 equals sTop1
    break;
  case 61 : isaID = 22; //22 equals sDownL
    break;
  case 62 : isaID = 23; //23 equals sStrangeL
    break;
  case 63 : isaID = 45; //25 equals sBottom2
    break;
  case 64 : isaID = 42; //42 equals sDownR
    break;
  case 65 : isaID = 43; //43 equals sStrangeR
    break;
  case 66 : isaID = 25; //43 equals sBottom1
    break;
  }
  if (antifl) isaID*=-1;
  return isaID;
}

Flavour Isajet_Fortran_Interface::IsaIDToFlavour(int isaID) {

  bool antifl  = 0;
  int sherpaID = 0;
  
  if (isaID<0) {
    isaID *= -1; 
    antifl = 1;
  }
  
  switch(isaID) {
  case 1   : sherpaID = 2; //up-quark
    break;
  case 2   : sherpaID = 1; //down-quark
    break;
  case 3   : sherpaID = 3; //strange-quark
    break;
  case 4   : sherpaID = 4; //charm-quark
    break;
  case 5   : sherpaID = 5; //bottom-quark
    break;
  case 6   : sherpaID = 6; //top-quark
    break;
  case 11  : sherpaID = 12; //nuelectron
    break;
  case 12  : sherpaID = 11; //electron
    break;
  case 13  : sherpaID = 14; //numu
    break;
  case 14  : sherpaID = 13; //muon
    break;
  case 15  : sherpaID = 16; //nutau
    break;
  case 16  : sherpaID = 15; //tau
    break;
  case 9   : sherpaID = 21; //gluon
    break;
  case 10  : sherpaID = 22; //photon
    break;
  case 80  : sherpaID = 23; //W-
    break;
  case 90  : sherpaID = 24; //Z0
    break;
  case 39  : sherpaID = 41; //Chargino1
    break;
  case 49  : sherpaID = 42; //Chargino2
    break;
  case 30  : sherpaID = 43; //Neutralino1
    break;
  case 40  : sherpaID = 44; //Neutralino2
    break;
  case 50  : sherpaID = 45; //Neutralino3
    break;
  case 60  : sherpaID = 46; //Neutralino4
    break;
  case 29  : sherpaID = 47; //gluino
    break;
  case 82  : sherpaID = 31; //h0
    break;  
  case 83  : sherpaID = 32; //H0
    break;
  case 84  : sherpaID = 33; //A0
    break;
  case 86  : sherpaID = 34; //Hmin
    break;
  case 31  : sherpaID = 81; //sNu1
     break;
  case 33  : sherpaID = 82; //sNu2
    break;
  case 35  : sherpaID = 83; //sNu3
    break;
  case 32  : sherpaID = 71; //sElectronL
    break;
  case 34  : sherpaID = 72; //sMuonL
    break;  
  case 56  : sherpaID = 73; //sTau2
    break;  
  case 52  : sherpaID = 74; //sElectronR
    break;  
  case 54  : sherpaID = 75; //sMuonR
    break;  
  case 36  : sherpaID = 76; //sTau1
    break;
  case 21  : sherpaID = 51; //sUpL
    break;
  case 24  : sherpaID = 52; //sCharmL
    break;
  case 46  : sherpaID = 53; //sTop2
    break;
  case 41  : sherpaID = 54; //sUpR
    break;
  case 44  : sherpaID = 55; //sCharmR
    break;
  case 26  : sherpaID = 56; //sTop1
    break;
  case 22  : sherpaID = 61; //sDownL 
    break;
  case 23  : sherpaID = 62; //sStrangeL 
    break;
  case 45  : sherpaID = 63; //sBottom2 
    break;
  case 42  : sherpaID = 64; //sDownR 
    break;
  case 43  : sherpaID = 65; //sStrangeR 
    break;
  case 25  : sherpaID = 66; //sBottom1 
    break;
  }
  Flavour flav = Flavour((kf_code)(sherpaID));
  if (sherpaID==23 || sherpaID==34) flav = flav.Bar(); 
  if (antifl) flav = flav.Bar();
  return flav;
}


void Isajet_Fortran_Interface::Decays(Flavour flav) {  
  int   mother,ndecays, * daughters, * ndaughters;
  float totalwidth, * partialwidths;
  daughters     = new int[300];
  ndaughters    = new int[100];
  partialwidths = new float[100];
 
  mother = FlavourToIsaID(flav);
  /*
    Take care for decays Hmin
    there the convention of Sherpa and Isasusy is opposite
    concerning the charge of the particles,
    same holds true for W-boson
  */    
  filldecays_(mother,totalwidth,ndecays,daughters,ndaughters,partialwidths);

  Decay_Table * dt = new Decay_Table(flav);
  m_widths.push_back(totalwidth);
  m_masses.push_back(flav.Mass());
  m_particles.push_back(flav);
  m_decays.push_back(dt);
  Decay_Channel * decay;

  for (int i=0;i<ndecays;i++) {
    decay = new Decay_Channel(flav);
    for (int j=0;j<3;j++) {
      if (daughters[3*i+j]!=0) decay->AddDecayProduct(IsaIDToFlavour(daughters[3*i+j]));
    }
    decay->SetWidth(partialwidths[i]);
    dt->AddDecayChannel(decay);
  }
  dt->SetWidthGenerator("Isajet");
  dt->Output();  
  flav.SetWidth(totalwidth);
  
  delete [] daughters;
  delete [] ndaughters;
  delete [] partialwidths;
}


void Isajet_Fortran_Interface::CharginoMasses()
{
  float MChi1,MChi2,gammaL,gammaR,ThX,ThY;
  CMatrix Zplus  = CMatrix(2);
  CMatrix Zminus = CMatrix(2);

  chargino_(MChi1,MChi2,gammaL,gammaR,ThX,ThY);

  MChi1 = -MChi1;
  MChi2 = -MChi2;

  int theta1=0,theta2=0;
  Flavour flav;
  flav = Flavour(kf_Chargino1);flav.SetMass(dabs(MChi1));
  if (MChi1<0) { flav.SetMassSign(-1); theta1 = 1; }
  flav = Flavour(kf_Chargino2);flav.SetMass(dabs(MChi2));
  if (MChi2<0) { flav.SetMassSign(-1); theta2 = 1; }

  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"Charginomasses :"<<std::endl
		<<Flavour(kf_Chargino1).MassSign()<<" * "<<Flavour(kf_Chargino1).Mass()
		<<", "
		<<Flavour(kf_Chargino2).MassSign()<<" * "<<Flavour(kf_Chargino2).Mass()<<std::endl;
  
  Zminus[0][0] = -(::sin(gammaL));
  Zminus[1][0] = -cos(gammaL);
  Zminus[0][1] = -ThX*cos(gammaL);
  Zminus[1][1] = ThX*(::sin(gammaL));
  
  Zplus[0][0] = -pow(-1.,theta1)*(::sin(gammaR));
  Zplus[1][0] = -pow(-1.,theta2)*cos(gammaR);
  Zplus[0][1] = -ThY*pow(-1.,theta1)*cos(gammaR);
  Zplus[1][1] = ThY*pow(-1.,theta2)*(::sin(gammaR));  
  
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z^+"),Zplus));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z^-"),Zminus));
}


void Isajet_Fortran_Interface::NeutralinoMasses()
{
  float * M   = new float[4];
  float * Mix = new float[16];

  
  neutralino_(M,Mix);
  for (short int i=0;i<4;i++) M[i] = -M[i];
  Matrix<4> M_Mix;
  for (short int i=0;i<4;i++) {
    for (short int j=0;j<4;j++) M_Mix[i][j] = Mix[i+j*4];
  }

  Flavour flav;
  flav = Flavour(kf_Neutralino1);flav.SetMass(dabs(M[0]));
  if (M[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_Neutralino2);flav.SetMass(dabs(M[1]));
  if (M[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_Neutralino3);flav.SetMass(dabs(M[2]));
  if (M[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_Neutralino4);flav.SetMass(dabs(M[3]));
  if (M[3]<0) flav.SetMassSign(-1);

  delete[] M;
  delete[] Mix;

  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"Neutralinomasses :"<<std::endl
		<<Flavour(kf_Neutralino1).MassSign()<<" * "<<Flavour(kf_Neutralino1).Mass()
		<<", "
		<<Flavour(kf_Neutralino2).MassSign()<<" * "<<Flavour(kf_Neutralino2).Mass()
		<<", "
		<<Flavour(kf_Neutralino3).MassSign()<<" * "<<Flavour(kf_Neutralino3).Mass()
		<<", "
		<<Flavour(kf_Neutralino4).MassSign()<<" * "<<Flavour(kf_Neutralino4).Mass()<<endl;

  CMatrix ZN = CMatrix(4);
  for (short int i=0;i<4;i++) {
    //transposed with respect to Richardson !
    ZN[0][i] = Complex(M_Mix[3][i],0.);
    ZN[1][i] = Complex(M_Mix[2][i],0.);
    ZN[2][i] = Complex(-M_Mix[1][i],0.);
    ZN[3][i] = Complex(-M_Mix[0][i],0.);        
  }
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z^N"),ZN));
}

void Isajet_Fortran_Interface::sNeutrinoMasses()
{
  float mn1,mn2,mn3;
  sneutrino_(mn1,mn2,mn3);
  
  Flavour flav;    
  flav = Flavour(kf_sNu1);flav.SetMass(dabs(mn1));
  if (mn1<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sNu2);flav.SetMass(dabs(mn2));
  if (mn2<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sNu3);flav.SetMass(dabs(mn3));
  if (mn3<0) flav.SetMassSign(-1);
  
  CMatrix ZNue = CMatrix(3);
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ZNue[i][j] = Complex(0.,0.);
    ZNue[i][i] = Complex(1.,0.);
  }    

  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl;
  msg_Tracking()<<"sNeutrinomasses :"<<std::endl
		<<Flavour(kf_sNu1).Mass()<<", "
		<<Flavour(kf_sNu2).Mass()<<", "
		<<Flavour(kf_sNu3).Mass()<<std::endl;

  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z_nu"),ZNue));
}

void Isajet_Fortran_Interface::HiggsMasses()
{
  float mh,mH,mA,mHC,tanb,alpha,vev;
  
  higgses_(mh,mH,mA,mHC,tanb,alpha,vev);

  // _v = vev;
  
  double cosb = sqrt(1./(1.+sqr(tanb)));
  double sinb = cosb*tanb;  
  
  CMatrix ZH  = CMatrix(2);
  ZH[0][0]    = Complex(sinb,0.);
  ZH[0][1]    = Complex(-cosb,0.);
  ZH[1][0]    = Complex(cosb,0.);
  ZH[1][1]    = Complex(sinb,0.);
    
  double sina = ::sin(alpha);
  double cosa = cos(alpha);

  CMatrix ZR  = CMatrix(2);
  ZR[0][0]    = Complex(cosa,0.);
  ZR[0][1]    = Complex(-sina,0.);
  ZR[1][0]    = Complex(sina,0.);
  ZR[1][1]    = Complex(cosa,0.);
  Flavour flav;
  flav = Flavour(kf_h0);flav.SetMass(mh);
  flav = Flavour(kf_H0);flav.SetMass(mH);
  flav = Flavour(kf_A0);flav.SetMass(mA);
  flav = Flavour(kf_Hmin);flav.SetMass(mHC);

  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"Higgs Masses :"<<std::endl
		<<"m_h0 = "<<Flavour(kf_h0).Mass()<<", m_H0 = "<<Flavour(kf_H0).Mass()
		<<"m_A0 = "<<Flavour(kf_A0).Mass()<<", m_H- = "<<Flavour(kf_Hmin).Mass()<<std::endl;

  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z_R"),ZR));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z_H"),ZH));
}

void Isajet_Fortran_Interface::sUpMasses()
{
  float * msups = new float[6];
  float thetaf,mu_1,_us22;
  sups_(msups,thetaf,mu_1,_us22);

  double mu    = mu_1;  
  double us22  = _us22;
  
  p_model->GetScalarConstants()->insert(std::make_pair(std::string("mu"),mu));

  double sinth = ::sin(thetaf);
  double costh = cos(thetaf);

  Flavour flav;    

  flav = Flavour(kf_sUpL);flav.SetMass(dabs(msups[0]));
  if (msups[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sCharmL);flav.SetMass(dabs(msups[1]));
  if (msups[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sTop2);flav.SetMass(dabs(msups[2]));
  if (msups[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sUpR);flav.SetMass(dabs(msups[3]));
  if (msups[3]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sCharmR);flav.SetMass(dabs(msups[4]));
  if (msups[4]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sTop1);flav.SetMass(dabs(msups[5]));
  if (msups[5]<0) flav.SetMassSign(-1);
  
  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"sQuark masses :"<<std::endl;
  for (short int i=0;i<5;++i) {
    msg_Tracking()<<"sUpquarks["<<i<<"] : "<<dabs(msups[i])<<" ,";
    if (i==2) msg_Tracking()<<endl;
  }
  msg_Tracking()<<"sUpquarks["<<5<<"] : "<<dabs(msups[5])<<" ,"<<std::endl;


  CMatrix Zu = CMatrix(6);
  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zu[i][j] = Complex(0.,0.);
    Zu[i][i] = Complex(1.,0.);
  }    
  Zu[2][2] = Complex(sinth,0.);
  Zu[2][5] = Complex(costh,0.);
  Zu[5][5] = Complex(sinth,0.);
  Zu[5][2] = Complex(-costh,0.);
  
  delete[] msups;

  CMatrix ws = CMatrix(3);
  CMatrix us = CMatrix(3);
  for (short int i=0;i<2;i++) {
    for (short int j=0;j<2;j++) ws[i][j] = Complex(0.,0.);   
  }    
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) us[i][j] = Complex(0.,0.);
    us[2][2] = Complex(us22,0.);  
  }    
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("ws"),ws));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("us"),us));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z_u"),Zu));
}

void Isajet_Fortran_Interface::sDownMasses()
{
  float* msdowns = new float[6];
  float thetaf,_ds22;
  sdowns_(msdowns,thetaf,_ds22);

  double ds22  = _ds22;
  double sinth = ::sin(thetaf);
  double costh = cos(thetaf);

  Flavour flav;    

  flav = Flavour(kf_sDownL);flav.SetMass(dabs(msdowns[0]));
  if (msdowns[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sStrangeL);flav.SetMass(dabs(msdowns[1]));
  if (msdowns[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sBottom2);flav.SetMass(dabs(msdowns[2]));
  if (msdowns[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sDownR);flav.SetMass(dabs(msdowns[3]));
  if (msdowns[3]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sStrangeR);flav.SetMass(dabs(msdowns[4]));
  if (msdowns[4]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sBottom1);flav.SetMass(dabs(msdowns[5]));
  if (msdowns[5]<0) flav.SetMassSign(-1);
  
  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"sQuark masses :"<<std::endl;
  for (short int i=0;i<5;++i) {
    msg_Tracking()<<"sDownquarks["<<i<<"] : "<<dabs(msdowns[i])<<" ,";
    if (i==2) msg_Tracking()<<endl;
  }
  msg_Tracking()<<"sDownquarks["<<5<<"] : "<<dabs(msdowns[5])<<" ,"<<std::endl;

  CMatrix Zd = CMatrix(6);
  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zd[i][j] = Complex(0.,0.);
    Zd[i][i] = Complex(1.,0.);
  }    
  Zd[2][2] = Complex(-sinth,0.);
  Zd[2][5] = Complex(costh,0.);
  Zd[5][5] = Complex(sinth,0.);
  Zd[5][2] = Complex(costh,0.);
  
  delete[] msdowns;

  CMatrix es = CMatrix(3);
  CMatrix ds = CMatrix(3);
  for (short int i=0;i<2;i++) {
    for (short int j=0;j<2;j++) es[i][j] = Complex(0.,0.);   
  }    
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ds[i][j] = Complex(0.,0.);
    ds[2][2] = Complex(ds22,0.);  
  }    
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("es"),es));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("ds"),ds));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z_d"),Zd));
}

void Isajet_Fortran_Interface::sLeptonMasses()
{
  float* msleptons = new float[6];
  float thetaf,mu_1,_ls22;
  sleptons_(msleptons,thetaf,mu_1,_ls22);

  double sinth = ::sin(thetaf);
  double costh = cos(thetaf);

  Flavour flav;    

  flav = Flavour(kf_sElectronL);flav.SetMass(dabs(msleptons[0]));
  if (msleptons[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sMuL);flav.SetMass(dabs(msleptons[1]));
  if (msleptons[1]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sTau2);flav.SetMass(dabs(msleptons[2]));
  if (msleptons[2]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sElectronR);flav.SetMass(dabs(msleptons[3]));
  if (msleptons[3]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sMuR);flav.SetMass(dabs(msleptons[4]));
  if (msleptons[4]<0) flav.SetMassSign(-1);
  flav = Flavour(kf_sTau1);flav.SetMass(dabs(msleptons[5]));
  if (msleptons[5]<0) flav.SetMassSign(-1);
  
  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"sLepton Masses :"<<std::endl;
  for (short int i=0;i<5;++i) {
    msg_Tracking()<<"sLeptons["<<i<<"] : "<<dabs(msleptons[i])<<" ,";
    if (i==2) msg_Tracking()<<endl;
  }
  msg_Tracking()<<"sLeptons["<<5<<"] : "<<dabs(msleptons[5])<<" ,"<<std::endl;
   
  CMatrix Zl = CMatrix(6);
  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zl[i][j] = Complex(0.,0.);
    Zl[i][i] = Complex(1.,0.);
  }    
  Zl[2][2] = Complex(sinth,0.);
  Zl[2][5] = Complex(costh,0.);
  Zl[5][5] = Complex(sinth,0.);
  Zl[5][2] = Complex(-costh,0.);
  

  for (short int i=0;i<6;i++) {
    for (short int j=0;j<6;j++) Zl[i][j] = 0.;
    Zl[i][i] = 1.;}    

  Zl[2][2] = sinth;
  Zl[2][5] = costh;
  Zl[5][5] = sinth;
  Zl[5][2] = -costh;
  
  delete[] msleptons;

  CMatrix ls = CMatrix(3);
  CMatrix ks = CMatrix(3);
  for (short int i=0;i<2;i++) {
    for (short int j=0;j<2;j++) ks[i][j] = Complex(0.,0.);   
  }    
  for (short int i=0;i<3;i++) {
    for (short int j=0;j<3;j++) ls[i][j] = Complex(0.,0.);
    ls[2][2] = Complex(double(_ls22),0.);  
  }    
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("ls"),ls));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("ks"),ks));
  p_model->GetComplexMatrices()->insert(std::make_pair(std::string("Z_l"),Zl));
}


void Isajet_Fortran_Interface::GluinoMasses()
{
  float _mGL;
  gluino_(_mGL);

  double mGL  = _mGL;
  
  Flavour flav;    

  flav = Flavour(kf_Gluino);flav.SetMass(dabs(mGL));
  if (mGL<0) flav.SetMassSign(-1);
  
  msg_Tracking()<<"--------------------------------------------------------------"<<std::endl
		<<"Gluino mass : "<<dabs(mGL)<<std::endl;
}
