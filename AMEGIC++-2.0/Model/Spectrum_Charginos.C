#include "Spectrum_Charginos.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Running_AlphaQED.H"
#include "MathTools.H"

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace std;

void Spectrum_Charginos::Interface(Isajet* isa)
{
  msg.Tracking()<<"=====================ISAJET==========================="<<endl;

  isa->Chargino(_Zplus,_Zminus);

  msg.Tracking()<<"Zplus: "<<endl;
  if (rpa.gen.Tracking()) _Zplus.MatrixOut();
  msg.Tracking()<<"Zminus: "<<endl;
  if (rpa.gen.Tracking()) _Zminus.MatrixOut();
  msg.Tracking()<<"======================================================"<<endl;
}

void Spectrum_Charginos::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  v     = dr.GetValue<double>("v");
  tanb  = dr.GetValue<double>("tan(beta)");
  m_2   = dr.GetValue<double>("m_2");
  mu    = dr.GetValue<double>("mu");

  CplEW.Init();
  Masses_LO();
}

void Spectrum_Charginos::Masses_LO()
{
  //   double e  = sqrt(4.*M_PI*(*aqed)(sqr(rpa.gen.Ecms())));
  double e  = sqrt(4.*M_PI*(aqed->AqedFixed()));

  double v1 = v*sqrt(1./(1.+sqr(tanb)));
  double v2 = v1*tanb;  

  Matrix<2> M;

  M[0][0] = 2.*m_2;
  M[0][1] = e*v2/sqrt(2.)/CplEW.SinThetaW();
  M[1][0] = e*v1/sqrt(2.)/CplEW.SinThetaW();
  M[1][1] = mu;
 

  double evalues[2];
  
  Matrix<2> MDM = M.Dagger()*M;
  MDM.DiagonalizeSort(evalues,_Zplus);

  Matrix<2> MMD = M*M.Dagger();
  MMD.DiagonalizeSort(evalues,_Zminus);

  /*
  Matrix<2> Dummy;

  Dummy[1][0] = 1.;
  Dummy[0][1] = 1.;

  _Zplus  = _Zplus*Dummy;
  _Zminus = _Zminus*Dummy;
  */


  Matrix<2> Dummy;

  Dummy = _Zminus.Dagger()*M*_Zplus;
  
  evalues[0] = Dummy[0][0];
  evalues[1] = Dummy[1][1];

  msg.Tracking()<<"Zplus: "<<endl;
  _Zplus.MatrixOut();
  msg.Tracking()<<"Zminus: "<<endl;
  _Zminus.MatrixOut();
  
  Flavour flav;
  flav = Flavour(kf::Chargino1);flav.SetMass(dabs(evalues[0]));
  if (evalues[0]<0) flav.SetMassSign(-1);
  flav = Flavour(kf::Chargino2);flav.SetMass(dabs(evalues[1]));
  if (evalues[1]<0) flav.SetMassSign(-1);
  
  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"Charginomasses :"<<endl;
  msg.Tracking()<<"m_Chi_1 = "<<Flavour(kf::Chargino1).Mass();
  msg.Tracking()<<" Sign: "<<Flavour(kf::Chargino1).MassSign();
  msg.Tracking()<<", m_Chi_2 = "<<Flavour(kf::Chargino2).Mass();
  msg.Tracking()<<" Sign: "<<Flavour(kf::Chargino2).MassSign()<<endl;
}

/*

This seems to be old ....

void Spectrum_Charginos::Isajet()
{
  string modelfilename=rpa.Get_Path()+string("/")+string("Isajet.dat");

  int parnum = 5;
  Data_Pointer** param = new Data_Pointer*[parnum];
  
  param[0] = new Data<Switch::code>;

  for (short int i=1;i<parnum;i++)
    param[i] = new Data<double>;

  param[0]->Set_Name(string("Isajet"));
  param[1]->Set_Name(string("mChi1"));
  param[2]->Set_Name(string("mChi2"));
  param[3]->Set_Name(string("gammaL"));
  param[4]->Set_Name(string("gammaR"));

  rpa.Read(modelfilename,param,parnum);

  Switch::code switch_tag;

  if (param[0]->Get_Value(switch_tag)==Switch::Off) return;

  double double_tag;

  Flavour flav;
  flav = Flavour(kf::Chargino1);flav.SetMass(param[1]->Get_Value(double_tag));
  flav = Flavour(kf::Chargino2);flav.SetMass(param[2]->Get_Value(double_tag));


  double gammaL = param[3]->Get_Value(double_tag); 
  double gammaR = param[4]->Get_Value(double_tag); 

  for (short int i=0;i<parnum;i++)
    delete param[i];
  delete[] param;

  double alpha = atan(tanb);

//  double e  = sqrt(4.*M_PI*(*aqed)(sqr(rpa.gen.Ecms())));
  double e  = sqrt(4.*M_PI*(aqed->AqedFixed()));
  double MW = e/(2.*CplEW.SinThetaW())*v;

  double Xi = sqr(4.*sqr(m_2)-sqr(mu))+4.*sqr(MW)
    *(sqr(MW)*sqr(cos(2.*alpha))+4.*sqr(m_2)+sqr(mu)+4.*m_2*mu*sin(2.*alpha));
								     
  double xminus = ((4.*sqr(m_2)-sqr(mu)-2.*sqr(MW)*cos(2.*alpha))-Xi)/
    (2.*sqrt(2.)*MW*(mu*sin(alpha)+2.*m_2*cos(alpha)));

  double yminus = ((4.*sqr(m_2)-sqr(mu)+2.*sqr(MW)*cos(2.*alpha))-Xi)/
    (2.*sqrt(2.)*MW*(mu*cos(alpha)+2.*m_2*sin(alpha)));
  
  int ThetaX = xminus>0 ? 1 : -1;
  int ThetaY = yminus>0 ? 1 : -1;

  double mplus = ThetaX*ThetaY*(cos(gammaR)*(mu*cos(gammaL)-e*v*sin(gammaL)*sin(alpha))-
				sin(gammaR)*(e*v*cos(gammaL)*cos(alpha)-2.*m_2*sin(gammaL)));

  double mminus = sin(gammaR)*(mu*sin(gammaL)+e*v*cos(gammaL)*sin(alpha))+
                  cos(gammaR)*(e*v*sin(gammaL)*cos(alpha)+2.*m_2*cos(gammaL));
  
  int ThetaP = mplus>0 ? 0 : 1;
  int ThetaM = mminus>0 ? 0 : 1;

  _Zminus[0][0] = -sin(gammaL);
  _Zminus[0][1] = -cos(gammaL);
  _Zminus[1][0] = -ThetaX*cos(gammaL);
  _Zminus[1][1] = ThetaX*sin(gammaL);

  _Zplus[0][0] = -sin(gammaR)*pow(-1,ThetaM);
  _Zplus[0][1] = -cos(gammaR)*pow(-1,ThetaP);
  _Zplus[1][0] = -ThetaY*cos(gammaR)*pow(-1,ThetaM);
  _Zplus[1][1] = ThetaY*sin(gammaR)*pow(-1,ThetaP);  

  _Zminus.MatrixOut();
  _Zplus.MatrixOut();
  msg.Tracking()<<"--------------------------------------------------------------"<<endl;
  msg.Tracking()<<"Charginomasses (IsaJet) :"<<endl;
  msg.Tracking()<<"m_Chi_1 = "<<Flavour(kf::Chargino1).Mass();
  msg.Tracking()<<", m_Chi_2 = "<<Flavour(kf::Chargino2).Mass()<<endl;
}



*/









