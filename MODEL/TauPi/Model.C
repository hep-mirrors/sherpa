#include "MODEL/Main/Model_Base.H"

namespace MODEL {

  class Standard_Model_TauPi: public Model_Base {
  private:

    int  m_ckmorder, m_dec_g4;

    void FixEWParameters();
    void FixCKM();

    void ParticleInit();

    void InitQEDVertices();
    void InitQCDVertices();
    void InitEWVertices();

  public :

    Standard_Model_TauPi(std::string,std::string);
    bool ModelInit(const PDF::ISR_Handler_Map& isr);
    void InitVertices();

  };

}

#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Standard_Model_TauPi,"TauPi",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,Standard_Model_TauPi>::
operator()(const Model_Arguments &args) const
{
  return new Standard_Model_TauPi(args.m_path,args.m_file);
}

void Getter<Model_Base,Model_Arguments,Standard_Model_TauPi>::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"The Standard Model with \\tau\\pi decays\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (values 0,1,3, EW input schemes, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ) (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- ALPHAQED_DEFAULT_SCALE (scale for alpha_QED default)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- CKMORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CABIBBO (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- A (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- RHO (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- ETA (Wolfenstein Eta)\n"
     <<setw(width+4)<<" "<<"}";
  str<<"Infrared continuation of alphaS:\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"- AS_FORM (values 0,1,2,3,10, see documentation)\n"
     <<setw(width+7)<<" "<<"- Q2_AS (corresponding infrared parameter, see documentation)\n"
     <<setw(width+4)<<" "<<"}";
}

Standard_Model_TauPi::Standard_Model_TauPi(string _dir,string _file) :
  Model_Base(_dir,_file,1)
{
  m_name="TauPi";
  ParticleInit();
  CustomContainerInit();
}

void Standard_Model_TauPi::ParticleInit()
{
  s_kftable[kf_none] = new ATOOLS::Particle_Info(kf_none,-1,0,0,0,0,-1,0,1,0,"no_particle","no_particle","no_particle", "no_particle", 1,1);
  //add SM particles
  //kf_code,mass,width,charge,strong,spin,majorana,take,stable,massive,idname,antiname,texname,antitexname
  s_kftable[kf_d]      = new Particle_Info(kf_d,0.01,.0,-1,3,1,0,1,1,0,"d","db", "d", "\\bar{d}");
  s_kftable[kf_u]      = new Particle_Info(kf_u,0.005,.0,2,3,1,0,1,1,0,"u","ub", "u", "\\bar{u}");
  s_kftable[kf_s]      = new Particle_Info(kf_s,0.2,.0,-1,3,1,0,1,1,0,"s","sb", "s", "\\bar{s}");
  s_kftable[kf_c]      = new Particle_Info(kf_c,1.42,.0,2,3,1,0,1,1,0,"c","cb", "c", "\\bar{c}");
  s_kftable[kf_b]      = new Particle_Info(kf_b,4.8,.0,-1,3,1,0,1,1,0,"b","bb", "b", "\\bar{b}");
  s_kftable[kf_t]      = new Particle_Info(kf_t,173.21,2.0,2,3,1,0,1,0,1,"t","tb", "t", "\\bar{t}");
  s_kftable[kf_e]      = new Particle_Info(kf_e,0.000511,.0,-3,0,1,0,1,1,0,"e-","e+", "e^{-}", "e^{+}");
  s_kftable[kf_nue]    = new Particle_Info(kf_nue,.0,.0,0,0,1,0,1,1,0,"ve","veb", "\\nu_{e}", "\\bar{\\nu}_{e}");
  s_kftable[kf_mu]     = new Particle_Info(kf_mu,.105,.0,-3,0,1,0,1,1,0,"mu-","mu+", "\\mu^{-}", "\\mu^{+}");
  s_kftable[kf_numu]   = new Particle_Info(kf_numu,.0,.0,0,0,1,0,1,1,0,"vmu","vmub", "\\nu_{\\mu}", "\\bar{\\nu}_{\\mu}");
  s_kftable[kf_tau]    = new Particle_Info(kf_tau,1.777,2.26735e-12,-3,0,1,0,1,0,0,"tau-","tau+", "\\tau^{-}", "\\tau^{+}");
  s_kftable[kf_nutau]  = new Particle_Info(kf_nutau,.0,.0,0,0,1,0,1,1,0,"vtau","vtaub", "\\nu_{\\tau}", "\\bar{\\nu}_{\\tau}");
  s_kftable[kf_gluon]  = new Particle_Info(kf_gluon,.0,.0,0,8,2,-1,1,1,0,"G","G", "G", "G");
  s_kftable[kf_photon] = new Particle_Info(kf_photon,.0,.0,0,0,2,-1,1,1,0,"P","P","\\gamma","\\gamma");
  s_kftable[kf_Z]      = new Particle_Info(kf_Z,91.1876,2.4952,0,0,2,-1,1,0,1,"Z","Z","Z","Z");
  s_kftable[kf_Wplus]  = new Particle_Info(kf_Wplus,80.385,2.085,3,0,2,0,1,0,1,"W+","W-","W^{+}","W^{-}");
  s_kftable[kf_h0]     = new Particle_Info(kf_h0,125.,0.00407,0,0,0,-1,1,0,1,"h0","h0","h_{0}","h_{0}");
  s_kftable[kf_gluon_qgc] = new Particle_Info(kf_gluon_qgc,0.0,0.0,0,8,4,-1,1,1,0,"G4","G4","G_{4}","G_{4}",1);
  s_kftable[kf_pi_plus] = new Particle_Info(kf_pi_plus,0.13957,2.5242e-17,3,0,1,1,"pi+","pi-");
  ReadParticleData();
  AddStandardContainers();
}

bool Standard_Model_TauPi::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  FixEWParameters();  
  FixCKM();
  SetAlphaQCD(isr,p_dataread->GetValue<double>("ALPHAS(MZ)",0.118));
  SetRunningFermionMasses();
  p_constants->insert(make_pair(string("F_PI"),
                                p_dataread->GetValue<double>("F_PI",0.0924)));
  
  ATOOLS::OutputParticles(msg->Info());
  ATOOLS::OutputContainers(msg->Info());
  OutputCKM();

  for (MODEL::ScalarNumbersMap::iterator it=p_numbers->begin();
       it!=p_numbers->end(); ++it) {
    DEBUG_INFO(it->first+" = "<<it->second);
  }
  for (MODEL::ScalarConstantsMap::iterator it=p_constants->begin();
       it!=p_constants->end(); ++it) {
    DEBUG_INFO(it->first+" = "<<it->second);
  }
  for (MODEL::ComplexConstantsMap::iterator it=p_complexconstants->begin();
       it!=p_complexconstants->end(); ++it) {
    DEBUG_INFO(it->first+" = "<<it->second);
  }


  return true;
}

void Standard_Model_TauPi::FixEWParameters()
{
  double alphaQED0;
  Complex csin2thetaW, ccos2thetaW, cvev, I(0.,1.);
  string yukscheme=p_dataread->GetValue<string>("YUKAWA_MASSES","Running");
  p_numbers->insert(make_pair(string("YukawaScheme"), yukscheme=="Running"));
  string widthscheme=p_dataread->GetValue<string>("WIDTH_SCHEME","CMS");
  p_numbers->insert(make_pair(string("WidthScheme"), widthscheme=="CMS"));
  int ewscheme=p_dataread->GetValue<int>("EW_SCHEME",1);
  double MW=Flavour(kf_Wplus).Mass(), GW=Flavour(kf_Wplus).Width();
  double MZ=Flavour(kf_Z).Mass(), GZ=Flavour(kf_Z).Width();
  double MH=Flavour(kf_h0).Mass(), GH=Flavour(kf_h0).Width();
  switch (ewscheme) {
  case 0:
    // all SM parameters given explicitly
    SetAlphaQEDByScale(p_dataread->GetValue<double>("ALPHAQED_DEFAULT_SCALE",sqr(MZ)));
    csin2thetaW=p_dataread->GetValue<double>("SIN2THETAW",0.23);
    ccos2thetaW=1.-csin2thetaW;
    cvev=p_dataread->GetValue<double>("VEV",246.);
    if (widthscheme=="CMS") THROW(not_implemented,"No CMS in EW_SCHEME=0");
    break;
  case 1: {
    // SM parameters given by alphaQED0, M_W, M_Z, M_H
    SetAlphaQEDByScale(p_dataread->GetValue<double>("ALPHAQED_DEFAULT_SCALE",sqr(MZ)));
    ccos2thetaW=sqr(MW/MZ);
    csin2thetaW=1.-ccos2thetaW;
    cvev=2.*MW*sqrt(csin2thetaW/(4.*M_PI*aqed->Default()));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=2.*sqrt(muW2*csin2thetaW/(4.*M_PI*aqed->Default()));
    }
    break;
  }
  case 2: {
    // SM parameters given by alphaQED(mZ), M_W, M_Z, M_H
    SetAlphaQED(1./p_dataread->GetValue<double>("1/ALPHAQED(MZ)",128.802));
    ccos2thetaW=sqr(MW/MZ);
    csin2thetaW=1.-ccos2thetaW;
    cvev=2.*MW*sqrt(csin2thetaW/(4.*M_PI*aqed->Default()));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=2.*sqrt(muW2*csin2thetaW/(4.*M_PI*aqed->Default()));
    }
    break;
  }
  case 3: {
    //gmu scheme -- inputs GF, M_W, M_Z, M_H
    double GF=p_dataread->GetValue<double>("GF",1.16639e-5);
    csin2thetaW=1.-sqr(MW/MZ);
    ccos2thetaW=1.-csin2thetaW;
    cvev=1./(pow(2.,0.25)*sqrt(GF));
    SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MW)*std::abs(csin2thetaW));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=1./(pow(2.,0.25)*sqrt(GF));
      break;
    }
    break;
  }
  case 4: {
    // FeynRules scheme, inputs: alphaQED, GF, M_Z, M_H
    SetAlphaQED(1./p_dataread->GetValue<double>("1/ALPHAQED(0)",137.03599976));
    double GF=p_dataread->GetValue<double>("GF",1.16639e-5);
    MW = sqrt(sqr(MZ)/2. + sqrt(pow(MZ,4)/4. - (aqed->Default()*M_PI*sqr(MZ))/(GF*sqrt(2.))));
    Flavour(kf_Wplus).SetMass(MW);

    csin2thetaW=1.-sqr(MW/MZ);
    ccos2thetaW=1.-csin2thetaW;
    cvev=1./(pow(2.,0.25)*sqrt(GF));

    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ)), muH2(MH*(MH-I*GH));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=1./(pow(2.,0.25)*sqrt(GF));
      break;
    }
    break;
  }
  default:
    THROW(not_implemented, "Unknown EW_SCHEME="+ToString(ewscheme));
    break;
  }
  p_complexconstants->insert(make_pair(string("ccos2_thetaW"),ccos2thetaW));
  p_complexconstants->insert(make_pair(string("csin2_thetaW"),csin2thetaW));
  p_complexconstants->insert(make_pair(string("cvev"), cvev));
}

void Standard_Model_TauPi::FixCKM()
{
  CMatrix CKM(3);
  for (int i=0;i<3;i++) {
    for (int j=i;j<3;j++) CKM[i][j] = CKM[j][i] = Complex(0.,0.);
    CKM[i][i] = Complex(1.,0.);
  }
  double Cabibbo=0.0,A=.8,rho,eta;
  m_ckmorder     = p_dataread->GetValue<int>("CKMORDER",0);
  if (m_ckmorder>0) {
    Cabibbo    = p_dataread->GetValue<double>("CABIBBO",0.22537);
    CKM[0][0] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[1][1] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[0][1] += Cabibbo * Complex( 1.,0.);
    CKM[1][0] += Cabibbo * Complex(-1.,0.);
  }
  if (m_ckmorder>1) {
    A          = p_dataread->GetValue<double>("A",0.814);
    CKM[1][2] += A*sqr(Cabibbo)  * Complex( 1.,0.);
    CKM[2][1] += A*sqr(Cabibbo)  * Complex(-1.,0.);
  }
  if (m_ckmorder>2) {
    eta        = p_dataread->GetValue<double>("ETA",0.353);
    rho        = p_dataread->GetValue<double>("RHO",0.117);
    CKM[0][2] += A*pow(Cabibbo,3) * Complex(rho,-eta);
    CKM[2][0] += A*pow(Cabibbo,3) * Complex(1.-rho,-eta);
  }

  ReadExplicitCKM(CKM);

  p_constants->insert(make_pair("CKM_DIMENSION",3));
  for (size_t i(0);i<3;++i)
    for (size_t j(0);j<3;++j)
      p_complexconstants->insert
	(make_pair("CKM_"+ToString(i)+"_"+ToString(j),CKM[i][j]));
  for (size_t i(0);i<3;++i)
    for (size_t j(0);j<3;++j)
      p_complexconstants->insert
	(make_pair("L_CKM_"+ToString(i)+"_"+ToString(j),i==j?1.0:0.0));
}

void Standard_Model_TauPi::InitVertices()
{
  InitQEDVertices();
  InitQCDVertices();
  InitEWVertices();
  Kabbala GF("G_F",1.0/sqrt(2.0)/std::abs(sqr(ComplexConstant("cvev"))));
  Kabbala f_pi("F_PI",ScalarConstant("F_PI"));
  Kabbala I("i",Complex(0.,1.)), rt2("\\sqrt{2}",sqrt(2.));
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_tau).Bar());
  m_v.back().AddParticle(Flavour(kf_nutau));
  m_v.back().AddParticle(Flavour(kf_pi_plus).Bar());
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("TauPiL");
  m_v.back().cpl.push_back(-I*f_pi*GF/rt2);
  m_v.back().order[1]=2;

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_nutau).Bar());
  m_v.back().AddParticle(Flavour(kf_tau));
  m_v.back().AddParticle(Flavour(kf_pi_plus));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("TauPiL");
  m_v.back().cpl.push_back(-I*f_pi*GF/rt2);
  m_v.back().order[1]=2;
}

void Standard_Model_TauPi::InitQEDVertices()
{
  if (!Flavour(kf_photon).IsOn()) return;
  Kabbala g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED")));
  Kabbala cpl=g1*Kabbala("i",Complex(0.,1.));
  for (short int i=1;i<17;++i) {
    if (i==7) i=11;
    Flavour flav((kf_code)i);
    if (flav.IsOn() && flav.Charge()) {
      Kabbala Q("Q_{"+flav.TexName()+"}",flav.Charge());
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_photon));
      m_v.back().Color.push_back
	(i>6?Color_Function(cf::None):
	 Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().cpl.push_back(cpl*Q);
      m_v.back().order[1]=1;
    } 
  }
}

void Standard_Model_TauPi::InitQCDVertices()
{
  if (!Flavour(kf_gluon).IsOn()) return;
  m_dec_g4 = p_dataread->GetValue<int>("DECOMPOSE_4G_VERTEX",1);
  Kabbala g3("g_3",sqrt(4.*M_PI*ScalarConstant("alpha_S")));
  Kabbala cpl0=g3*Kabbala("i",Complex(0.,1.));
  for (short int i=1;i<=6;++i) {
    Flavour flav((kf_code)i);
    if (!flav.IsOn()) continue; 
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_gluon));
    m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
    m_v.back().Lorentz.push_back("FFV");
    m_v.back().cpl.push_back(cpl0);
    m_v.back().order[0]=1;
  }
  Kabbala cpl1=-g3;
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<3;++i) m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().Color.push_back(Color_Function(cf::F,1,2,3));
  m_v.back().Lorentz.push_back("VVV");
  m_v.back().cpl.push_back(cpl1);
  m_v.back().order[0]=1;
  if (m_dec_g4) {
    m_v.push_back(Single_Vertex());
    for (size_t i(0);i<2;++i) m_v.back().AddParticle(Flavour(kf_gluon));
    m_v.back().AddParticle(Flavour(kf_gluon_qgc));
    m_v.back().Color.push_back(Color_Function(cf::F,1,2,3));
    m_v.back().Lorentz.push_back("VVP");
    m_v.back().cpl.push_back(cpl1);
    m_v.back().order[0]=1;
  }
  Kabbala cpl2=g3*g3*Kabbala("i",Complex(0.,1.)); 
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<4;++i) m_v.back().AddParticle(Flavour(kf_gluon));
  for (size_t i(0);i<3;++i) m_v.back().cpl.push_back(cpl2);
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,2,new Color_Function(cf::F,3,4,-1)));
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,3,new Color_Function(cf::F,2,4,-1)));
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,4,new Color_Function(cf::F,2,3,-1)));
  m_v.back().Lorentz.push_back("VVVVA");
  m_v.back().Lorentz.push_back("VVVVB");
  m_v.back().Lorentz.push_back("VVVVC");
  m_v.back().order[0]=2;
  if (m_dec_g4) m_v.back().dec=-1;
}

void Standard_Model_TauPi::InitEWVertices()
{
  Kabbala two(Kabbala("2",2.0)), three(Kabbala("3",3.0));
  Kabbala I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0));
  Kabbala g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED")));
  Kabbala sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW")));
  Kabbala costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW")));
  Kabbala g2(g1/sintW), vev("v_{EW}",ComplexConstant("cvev"));
  if (Flavour(kf_Wplus).IsOn()) {
    Kabbala cpl=I/rt2*g2;
    for (short int i=1;i<17;i+=2) {
      if (i==7) i=11;
      Flavour flav1((kf_code)i);
      if (!flav1.IsOn()) continue;
      for (short int j=2;j<18;j+=2) {
	if (j==8) j=12;
	if ((i<10 && j>10) || (i>10 && j<10)) continue;
	Flavour flav2((kf_code)j);
	if (!flav2.IsOn()) continue;
	std::string ckmstr=(i<10?"CKM_":"L_CKM_")+
	  ToString(((i%10)-1)/2)+"_"+ToString((j%10)/2-1);
	Kabbala ckm(ckmstr,ComplexConstant(ckmstr));
	if (std::abs(ckm.Value())==0.0) continue;
	m_v.push_back(Single_Vertex());
	m_v.back().AddParticle(flav1.Bar());
	m_v.back().AddParticle(flav2);
	m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
	m_v.back().Color.push_back
	  (i>6?Color_Function(cf::None):
	   Color_Function(cf::D,1,2));
	m_v.back().Lorentz.push_back("FFVL");
	m_v.back().cpl.push_back(cpl*ckm);
	m_v.back().order[1]=1;
	m_v.push_back(Single_Vertex());
	m_v.back().AddParticle(flav2.Bar());
	m_v.back().AddParticle(flav1);
	m_v.back().AddParticle(Flavour(kf_Wplus));
	m_v.back().Color.push_back
	  (i>6?Color_Function(cf::None):
	   Color_Function(cf::D,1,2));
	m_v.back().Lorentz.push_back("FFVL");
	m_v.back().cpl.push_back(cpl*ckm);
	m_v.back().order[1]=1;
      } 
    }
  }
  if (Flavour(kf_Z).IsOn()) {
    for (short int i=1;i<17;++i) {
      if (i==7) i=11;
      Flavour flav((kf_code)i);
      if (!flav.IsOn()) continue;
      Kabbala Q("Q_{"+flav.TexName()+"}",flav.Charge());
      Kabbala W("T_{"+flav.TexName()+"}",flav.IsoWeak());
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_Z));
      m_v.back().Color.push_back
	(i>6?Color_Function(cf::None):
	 Color_Function(cf::D,1,2));
      m_v.back().Color.push_back
	(i>6?Color_Function(cf::None):
	 Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().cpl.push_back(I/costW*(-Q*sintW+W/sintW)*g1);
      m_v.back().cpl.push_back(-I/costW*Q*sintW*g1);
      m_v.back().order[1]=1;
    } 
  }
  if (Flavour(kf_h0).IsOn()) {
    Kabbala cpl(-I/vev);
    for (short int i=1;i<17;++i) {
      if (i==7) i=11;
      Flavour flav((kf_code)i);
      if (!flav.IsOn() || flav.Yuk()==0.0) continue;
      double m=(ScalarNumber("YukawaScheme")==0)?flav.Yuk():
	ScalarFunction("m"+flav.IDName(),sqr(Flavour(kf_h0).Mass(true)));
      Kabbala M;
      if (ScalarNumber("WidthScheme")!=0)
        M=Kabbala("M_{"+flav.TexName()+"}(m_h^2)",
		  sqrt(m*m-Complex(0.0,m*flav.Width())));
      else M=Kabbala("M_{"+flav.TexName()+"}(m_h^2)",m);
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().Color.push_back
	(i>6?Color_Function(cf::None):
	 Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFS");
      m_v.back().cpl.push_back(cpl*M);
      m_v.back().order[1]=1;
    } 
  }
  if (Flavour(kf_Wplus).IsOn()) {
    if (Flavour(kf_photon).IsOn()) {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
      m_v.back().AddParticle(Flavour(kf_Wplus));
      m_v.back().AddParticle(Flavour(kf_photon));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVV");
      m_v.back().cpl.push_back(I*g1);
      m_v.back().order[1]=1;
    }      
    if (Flavour(kf_Z).IsOn()) {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
      m_v.back().AddParticle(Flavour(kf_Wplus));
      m_v.back().AddParticle(Flavour(kf_Z));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVV");
      m_v.back().cpl.push_back(I*g2*costW);
      m_v.back().order[1]=1;
    }
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
    m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
    m_v.back().AddParticle(Flavour(kf_Wplus));
    m_v.back().AddParticle(Flavour(kf_Wplus));
    m_v.back().cpl.push_back(-I*g2*g2);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVVV");
    m_v.back().order[1]=2;
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
    m_v.back().AddParticle(Flavour(kf_Wplus));
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().cpl.push_back(I*g1*g1);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVVV");
    m_v.back().order[1]=2;
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
    m_v.back().AddParticle(Flavour(kf_Wplus));
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(kf_Z));
    m_v.back().cpl.push_back(I*g1*g2*costW);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVVV");
    m_v.back().order[1]=2;
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
    m_v.back().AddParticle(Flavour(kf_Wplus));
    m_v.back().AddParticle(Flavour(kf_Z));
    m_v.back().AddParticle(Flavour(kf_Z));
    m_v.back().cpl.push_back(I*g2*g2*costW*costW);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVVV");
    m_v.back().order[1]=2;
  }
  if (Flavour(kf_h0).IsOn()) {
    if (Flavour(kf_Wplus).IsOn()) {
      Kabbala M("M_W",Flavour(kf_Wplus).Mass()), cpl;
      if (ScalarNumber("WidthScheme")!=0) {
	Kabbala G("\\Gamma_W",Flavour(kf_Wplus).Width());
	M=Kabbala("M_W",sqrt((M*M-I*G*M).Value()));
      }
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
      m_v.back().AddParticle(Flavour(kf_Wplus));
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().cpl.push_back(I*g2*M);
      m_v.back().order[1]=1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(Flavour(kf_Wplus).Bar());
      m_v.back().AddParticle(Flavour(kf_Wplus));
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().cpl.push_back(I*g2*g2/two);
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order[1]=2;
    }
    if (Flavour(kf_Z).IsOn()) {
      Kabbala M("M_Z",Flavour(kf_Z).Mass()), cpl;
      if (ScalarNumber("WidthScheme")!=0) {
	Kabbala G("\\Gamma_Z",Flavour(kf_Z).Width());
	M=Kabbala("M_Z",sqrt((M*M-I*G*M).Value()));
      }
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(Flavour(kf_Z));
      m_v.back().AddParticle(Flavour(kf_Z));
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().cpl.push_back(I*g2*M/costW);
      m_v.back().order[1]=1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(Flavour(kf_Z));
      m_v.back().AddParticle(Flavour(kf_Z));
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().AddParticle(Flavour(kf_h0));
      m_v.back().cpl.push_back(I*g2*g2/(costW*costW*two));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order[1]=2;
    }
    Kabbala M("M_H",Flavour(kf_h0).Mass()), cpl;
    if (ScalarNumber("WidthScheme")!=0) {
      Kabbala G("\\Gamma_H",Flavour(kf_h0).Width());
      M=Kabbala("M_H",sqrt((M*M-I*G*M).Value()));
    }
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSS");
    m_v.back().cpl.push_back(-I*M*M*three/vev);
    m_v.back().order[1]=1;
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().AddParticle(Flavour(kf_h0));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSSS");
    m_v.back().cpl.push_back(-I*M*M*three/(vev*vev));
    m_v.back().order[1]=2;
  }
}
