#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace MODEL {

  class Standard_ModelGS: public Model_Base {
  private:

    int  m_ckmorder, m_dec_g4;

    // void RegisterDefaults() const override;
    
    void FixEWParameters();
    void FixCKM();

    void ParticleInit() override;

    void InitQEDVertices();
    void InitQCDVertices();
    void InitEWVertices();
    void InitGoldstoneVertices();

    // util functions for phiplus vertices
    void GplusFermions();
    void GplusPhotonsV();
    void GplusZsV();
    void GplusHiggs();
    void GplusW();
    void GplusGold();

    // util functions for chi vertices
    void G0Fermions();
    void G0PhotonsV();
    void G0ZsV();
    void G0Higgs();
    void G0W();
    void G0Gold();
    
  public :

    Standard_ModelGS();
    bool ModelInit(const PDF::ISR_Handler_Map& isr) override;
    void InitVertices() override;
    void ResetVerticesWithEWParameters(const EWParameters&) override;
    void ClearInteractionModel();
    size_t IndexOfOrderKey(const std::string& key) const override;

  };

}

#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Standard_ModelGS,"SMGold",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,Standard_ModelGS>::
operator()(const Model_Arguments &args) const
{
  return new Standard_ModelGS();
}

void Getter<Model_Base,Model_Arguments,Standard_ModelGS>::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"The Standard Model\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"# possible parameters in yaml configuration [usage: \"keyword: value\"]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (values 0,1,3, EW input schemes, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ) (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- ALPHAQED_DEFAULT_SCALE (scale for alpha_QED default)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- CKM_ORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CKM_CABIBBO (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- CKM_A (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- CKM_RHO (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- CKM_ETA (Wolfenstein Eta)\n"
     <<setw(width+7)<<" "<<"- CKM_ELEMENT[<i>][<j>] (explicit value for element, supersedes parametrisation)\n"
     <<setw(width+4)<<" "<<"}";
  str<<"Infrared continuation of alphaS:\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"- AS_FORM (values 0,1,2,3,10, see documentation)\n"
     <<setw(width+7)<<" "<<"- Q2_AS (corresponding infrared parameter, see documentation)\n"
     <<setw(width+4)<<" "<<"}";
}

Standard_ModelGS::Standard_ModelGS() :
  Model_Base(true)
{
  m_name="SMGold";
  ParticleInit();
  RegisterDefaults();
  AddStandardContainers();
  CustomContainerInit();
}

void Standard_ModelGS::ParticleInit()
{
  // TODO: port this to all other models, or make sure in some other (not
  // model-specific) way, that the creation of more than one model leads to a
  // lot of memory leaks
  if (s_kftable.find(kf_none) != s_kftable.end()) {
    return;
  }
  s_kftable[kf_none] = new ATOOLS::Particle_Info(kf_none,-1,0,0,0,0,-1,0,1,0,"no_particle","no_particle","no_particle", "no_particle", 1,1);
  //add SM particles
  //kf_code,mass,width,charge,strong,spin,majorana,take,stable,massive,idname,antiname,texname,antitexname
  s_kftable[kf_d]         = new Particle_Info(kf_d,0.01,.0,-1,3,1,0,1,1,0,"d","db", "d", "\\bar{d}");
  s_kftable[kf_u]         = new Particle_Info(kf_u,0.005,.0,2,3,1,0,1,1,0,"u","ub", "u", "\\bar{u}");
  s_kftable[kf_s]         = new Particle_Info(kf_s,0.2,.0,-1,3,1,0,1,1,0,"s","sb", "s", "\\bar{s}");
  s_kftable[kf_c]         = new Particle_Info(kf_c,1.42,.0,2,3,1,0,1,1,0,"c","cb", "c", "\\bar{c}");
  s_kftable[kf_b]         = new Particle_Info(kf_b,4.8,.0,-1,3,1,0,1,1,0,"b","bb", "b", "\\bar{b}");
  s_kftable[kf_t]         = new Particle_Info(kf_t,173.21,2.0,2,3,1,0,1,0,1,"t","tb", "t", "\\bar{t}");
  s_kftable[kf_e]         = new Particle_Info(kf_e,0.000511,.0,-3,0,1,0,1,1,0,"e-","e+", "e^{-}", "e^{+}");
  s_kftable[kf_nue]       = new Particle_Info(kf_nue,.0,.0,0,0,1,0,1,1,0,"ve","veb", "\\nu_{e}", "\\bar{\\nu}_{e}");
  s_kftable[kf_mu]        = new Particle_Info(kf_mu,.105,.0,-3,0,1,0,1,1,0,"mu-","mu+", "\\mu^{-}", "\\mu^{+}");
  s_kftable[kf_numu]      = new Particle_Info(kf_numu,.0,.0,0,0,1,0,1,1,0,"vmu","vmub", "\\nu_{\\mu}", "\\bar{\\nu}_{\\mu}");
  s_kftable[kf_tau]       = new Particle_Info(kf_tau,1.777,2.26735e-12,-3,0,1,0,1,0,0,"tau-","tau+", "\\tau^{-}", "\\tau^{+}");
  s_kftable[kf_nutau]     = new Particle_Info(kf_nutau,.0,.0,0,0,1,0,1,1,0,"vtau","vtaub", "\\nu_{\\tau}", "\\bar{\\nu}_{\\tau}");
  s_kftable[kf_gluon]     = new Particle_Info(kf_gluon,.0,.0,0,8,2,-1,1,1,0,"G","G", "G", "G");
  s_kftable[kf_photon]    = new Particle_Info(kf_photon,.0,.0,0,0,2,-1,1,1,0,"P","P","\\gamma","\\gamma");
  s_kftable[kf_Z]         = new Particle_Info(kf_Z,91.1876,2.4952,0,0,2,-1,1,0,1,"Z","Z","Z","Z");
  s_kftable[kf_Wplus]     = new Particle_Info(kf_Wplus,80.385,2.085,3,0,2,0,1,0,1,"W+","W-","W^{+}","W^{-}");
  s_kftable[kf_h0]        = new Particle_Info(kf_h0,125.,0.00407,0,0,0,-1,1,0,1,"h0","h0","h_{0}","h_{0}");
  s_kftable[kf_gluon_qgc] = new Particle_Info(kf_gluon_qgc,0.0,0.0,0,8,4,-1,1,1,0,"G4","G4","G_{4}","G_{4}",1);
  
  s_kftable[kf_phiplus]   = new Particle_Info(kf_phiplus,80.385,2.085,3,0,0,0,1,0,1,"phi+","phi-","\\phi^{+}","\\phi^{-}");
  s_kftable[kf_chi]       = new Particle_Info(kf_chi,91.1876,2.4952,0,0,0,-1,1,0,1,"chi","chi","\\chi","\\chi");
  
  ReadParticleData();
}

bool Standard_ModelGS::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  FixEWParameters();  
  FixCKM();
  Settings& s = Settings::GetMainSettings();
  SetAlphaQCD(isr, s["ALPHAS(MZ)"].Get<double>());
  SetRunningFermionMasses();
  SetRunningBosonMasses();
  ATOOLS::OutputParticles(msg->Info());
  ATOOLS::OutputContainers(msg->Info());
  OutputCKM();
  for (MODEL::ScalarNumbersMap::iterator it=p_numbers->begin();
       it!=p_numbers->end();++it) DEBUG_INFO(it->first+" = "<<it->second);
  for (MODEL::ScalarConstantsMap::iterator it=p_constants->begin();
       it!=p_constants->end();++it) DEBUG_INFO(it->first+" = "<<it->second);
  for (MODEL::ComplexConstantsMap::iterator it=p_complexconstants->begin();
       it!=p_complexconstants->end();++it) DEBUG_INFO(it->first+" = "<<it->second);
  return true;
}

void Standard_ModelGS::FixEWParameters()
{
  Settings& s = Settings::GetMainSettings();
  Complex csin2thetaW, ccos2thetaW, cvev, I(0.,1.);
  string yukscheme = s["YUKAWA_MASSES"].Get<string>();
  p_numbers->insert(make_pair(string("YukawaScheme"), yukscheme=="Running"));
  string widthscheme = s["WIDTH_SCHEME"].Get<string>();
  p_numbers->insert(make_pair(string("WidthScheme"), widthscheme=="CMS"));
  int ewscheme = s["EW_SCHEME"].Get<int>();
  int ewrenscheme = s["EW_REN_SCHEME"].Get<int>();
  double MW=Flavour(kf_Wplus).Mass(), GW=Flavour(kf_Wplus).Width();
  double MZ=Flavour(kf_Z).Mass(), GZ=Flavour(kf_Z).Width();
  double MH=Flavour(kf_h0).Mass(), GH=Flavour(kf_h0).Width();
  std::string ewschemename(""),ewrenschemename("");
  switch (ewscheme) {
  case 0:
    // all SM parameters given explicitly
    ewschemename="user-defined, input: all parameters";
    SetAlphaQEDByScale(s["ALPHAQED_DEFAULT_SCALE"].Get<double>());
    csin2thetaW = s["SIN2THETAW"].Get<double>();
    ccos2thetaW=1.-csin2thetaW;
    cvev = s["VEV"].Get<double>();
    break;
  case 1: {
    // SM parameters given by alphaQED0, M_W, M_Z, M_H
    ewschemename="alpha(0) scheme, input: 1/\\alphaQED(0), m_W, m_Z, m_h, widths";
    SetAlphaQEDByScale(s["ALPHAQED_DEFAULT_SCALE"].Get<double>());
    ccos2thetaW=sqr(MW/MZ);
    csin2thetaW=1.-ccos2thetaW;
    cvev=2.*MW*sqrt(csin2thetaW/(4.*M_PI*aqed->Default()));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=2.*sqrt(muW2*csin2thetaW/(4.*M_PI*aqed->Default()));
    }
    break;
  }
  case 2: {
    // SM parameters given by alphaQED(mZ), M_W, M_Z, M_H
    ewschemename="alpha(m_Z) scheme, input: 1/\\alphaQED(m_Z), m_W, m_Z, m_h, widths";
    SetAlphaQEDByInput("1/ALPHAQED(MZ)");
    ccos2thetaW=sqr(MW/MZ);
    csin2thetaW=1.-ccos2thetaW;
    cvev=2.*MW*sqrt(csin2thetaW/(4.*M_PI*aqed->Default()));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=2.*sqrt(muW2*csin2thetaW/(4.*M_PI*aqed->Default()));
    }
    break;
  }
  case 3: {
    // Gmu scheme
    ewschemename="Gmu scheme, input: GF, m_W, m_Z, m_h, widths";
    double GF = s["GF"].Get<double>();
    csin2thetaW=1.-sqr(MW/MZ);
    ccos2thetaW=1.-csin2thetaW;
    cvev=1./(pow(2.,0.25)*sqrt(GF));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(MZ*(MZ-I*GZ));
      ccos2thetaW=muW2/muZ2;
      csin2thetaW=1.-ccos2thetaW;
      cvev=1./(pow(2.,0.25)*sqrt(GF));
      const size_t aqedconv{ s["GMU_CMS_AQED_CONVENTION"].Get<size_t>() };
      switch (aqedconv) {
      case 0:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::abs(muW2*csin2thetaW));
        break;
      case 1:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muW2*csin2thetaW));
        break;
      case 2:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muW2)*std::real(csin2thetaW));
        break;
      case 3 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MW)*std::abs(csin2thetaW));
        break;
      case 4 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MW)*(1.-sqr(MW/MZ)));
        break;
      default:
        THROW(not_implemented,"\\alpha_QED convention not implemented.");
      }
    } else if (widthscheme=="Fixed") {
      if (csin2thetaW.imag()!=0.0) THROW(fatal_error,"sin^2(\\theta_w) not real.");
      SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MW)*std::abs(csin2thetaW));
    }
    break;
  }
  case 4: {
    // DY scheme
    ewschemename="DY scheme, input: 1/\\alphaQED(m_Z), sin^2(theta_W), m_Z, m_h, widths";
    SetAlphaQEDByInput("1/ALPHAQED(MZ)");
    csin2thetaW = s["SIN2THETAW"].Get<double>();
    ccos2thetaW=1.-csin2thetaW;
    MW=MZ*sqrt(ccos2thetaW.real());
    Flavour(kf_Wplus).SetMass(MW);
    cvev=2.*MZ*sqrt(ccos2thetaW*csin2thetaW/(4.*M_PI*aqed->Default()));

    if (widthscheme=="CMS") {
      // now also the W width is defined by the tree-level relations
      Complex muW2(0.,0.), muZ2(MZ*(MZ-I*GZ));
      muW2=muZ2*ccos2thetaW;
      MW=sqrt(muW2.real());
      GW=-muW2.imag()/MW;
      Flavour(kf_Wplus).SetMass(MW);
      Flavour(kf_Wplus).SetWidth(GW);
      cvev=2.*sqrt(muZ2*ccos2thetaW*csin2thetaW/(4.*M_PI*aqed->Default()));
      break;
    }
    break;
  }
  case 5: {
    // CDY scheme
    ewschemename="CDY scheme, input: 1/\\alphaQED(m_W), sin^2(theta_W), m_W, m_h, widths";
    SetAlphaQEDByInput("1/ALPHAQED(MW)");
    csin2thetaW = s["SIN2THETAW"].Get<double>();
    ccos2thetaW=1.-csin2thetaW;
    MZ=MW/sqrt(ccos2thetaW.real());
    Flavour(kf_Z).SetMass(MZ);
    cvev=2.*MW*sqrt(csin2thetaW/(4.*M_PI*aqed->Default()));

    if (widthscheme=="CMS") {
      // now also the W width is defined by the tree-level relations
      Complex muW2(MW*(MW-I*GW)), muZ2(0.,0.);
      muZ2=muW2/ccos2thetaW;
      MZ=sqrt(muZ2.real());
      GZ=-muZ2.imag()/MZ;
      Flavour(kf_Z).SetMass(MZ);
      Flavour(kf_Z).SetWidth(GZ);
      cvev=2.*sqrt(muW2*csin2thetaW/(4.*M_PI*aqed->Default()));
      break;
    }
    break;
  }
  case 10: {
    // FeynRules scheme, inputs: alphaQED, GF, M_Z, M_H
    ewschemename="FeynRules scheme, input: 1/\\alphaQED(0), GF, m_Z, m_h, widths";
    SetAlphaQED(1./s["1/ALPHAQED(0)"].Get<double>());
    double GF = s["GF"].Get<double>();
    MW=sqrt(sqr(MZ)/2.+sqrt(pow(MZ,4)/4.
                            -(aqed->Default()*M_PI*sqr(MZ))/(GF*sqrt(2.))));
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
  switch (ewrenscheme) {
  case 1:
    ewrenschemename="alpha(0)";
    break;
  case 2:
    ewrenschemename="alpha(m_Z)";
    break;
  case 3:
    ewrenschemename="alpha(Gmu)";
    break;
  default:
    msg_Info()<<"Unknown EW_REN_SCHEME="<<ewrenscheme<<", resetting to 3."
              <<std::endl;
    ewrenscheme=3;
    ewrenschemename="alpha(Gmu)";
    break;
  }

  msg_Info()<<METHOD<<"() {"<<std::endl;
  msg_Info()<<"  Input scheme: "<<ewscheme<<std::endl;
  msg_Info()<<"                "<<ewschemename<<std::endl;
  msg_Info()<<"  Ren. scheme:  "<<ewrenscheme<<std::endl;
  msg_Info()<<"                "<<ewrenschemename<<std::endl;
  msg_Info()<<"  Parameters:   sin^2(\\theta_W) = "<<csin2thetaW.real()
            <<(csin2thetaW.imag()!=0.?(csin2thetaW.imag()>0?" + ":" - ")
                                       +ToString(abs(csin2thetaW.imag()),
                                                 msg->Precision())+" i"
                                     :"")<<std::endl;
  msg_Info()<<"                vev             = "<<cvev.real()
            <<(cvev.imag()!=0.?(cvev.imag()>0?" + ":" - ")
                                       +ToString(abs(cvev.imag()),
                                                 msg->Precision())+" i"
                                     :"")<<std::endl;
  msg_Info()<<"}"<<std::endl;
  aqed->PrintSummary();
  p_complexconstants->insert(make_pair(string("ccos2_thetaW"),ccos2thetaW));
  p_complexconstants->insert(make_pair(string("csin2_thetaW"),csin2thetaW));
  p_complexconstants->insert(make_pair(string("cvev"), cvev));
  rpa->gen.SetVariable("EW_SCHEME",ToString(ewscheme));
  rpa->gen.SetVariable("EW_REN_SCHEME",ToString(ewrenscheme));
}

void Standard_ModelGS::FixCKM()
{
  auto s = Settings::GetMainSettings()["CKM"];
  CMatrix CKM(3);
  for (int i=0;i<3;i++) {
    for (int j=i;j<3;j++) CKM[i][j] = CKM[j][i] = Complex(0.,0.);
    CKM[i][i] = Complex(1.,0.);
  }
  double Cabibbo=0.0,A=.8,rho,eta;
  m_ckmorder     = s["Order"].Get<int>();
  if (m_ckmorder>0) {
    Cabibbo    = s["Cabibbo"].Get<double>();
    CKM[0][0] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[1][1] += sqr(Cabibbo)/2. * Complex(-1.,0.);
    CKM[0][1] += Cabibbo * Complex( 1.,0.);
    CKM[1][0] += Cabibbo * Complex(-1.,0.);
  }
  if (m_ckmorder>1) {
    A          = s["A"].Get<double>();
    CKM[1][2] += A*sqr(Cabibbo)  * Complex( 1.,0.);
    CKM[2][1] += A*sqr(Cabibbo)  * Complex(-1.,0.);
  }
  if (m_ckmorder>2) {
    eta        = s["Eta"].Get<double>();
    rho        = s["Rho"].Get<double>();
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

void Standard_ModelGS::InitVertices()
{
  InitQEDVertices();
  InitQCDVertices();
  InitEWVertices();
  InitGoldstoneVertices();
}

void Standard_ModelGS::InitQEDVertices()
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

void Standard_ModelGS::InitQCDVertices()
{
  Settings& s = Settings::GetMainSettings();
  if (!Flavour(kf_gluon).IsOn()) return;
  m_dec_g4 = s["DECOMPOSE_4G_VERTEX"].Get<int>();
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
    m_v.back().dec=1;
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

void Standard_ModelGS::InitEWVertices()
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
      double m=ScalarConstant("m"+flav.IDName());
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
      Kabbala M("M_W",ScalarConstant("mW+")), cpl;
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
      Kabbala M("M_Z",ScalarConstant("mZ")), cpl;
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
    Kabbala M("M_H",ScalarConstant("mh0")), cpl;
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

void Standard_ModelGS::InitGoldstoneVertices()
{
  DEBUG_FUNC(METHOD);
  // phi vertexes
  if (Flavour(kf_phiplus).IsOn()) {
    // all fermion vertices
    //   GplusFermions();
    if (Flavour(kf_photon).IsOn()) {
      // all photon vertices
      GplusPhotonsV();
    }
    if (Flavour(kf_Z).IsOn()) {
      // all z vertices + zgamma
      GplusZsV();
    }
    // phi+ phi- -> h0 (EQ 96)     
    if (Flavour(kf_h0).IsOn()) {
      GplusHiggs();
    } // end kf_h0
    // phi± h/chi -> W±  and phi± W± -> A/Z     
    if (Flavour(kf_Wplus).IsOn()) {
      GplusW();
    } // end kf_Wplus
    GplusGold();
  } // end Flavour(kf_phiplus).IsOn()
  // chi vertexes
  if (Flavour(kf_chi).IsOn()) {
    G0Fermions();
    if (Flavour(kf_photon).IsOn()) {
      G0PhotonsV();
    }
    if (Flavour(kf_Z).IsOn()) {
      G0ZsV();
    }
    if (Flavour(kf_h0).IsOn()) {
      G0Higgs();
    }
    if (Flavour(kf_Wplus).IsOn()) {
      G0W();
    }
    G0Gold();
  } // end Flavour(kf_chi).IsOn()
}

void Standard_ModelGS::GplusFermions()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala two(Kabbala("2",2.0)), I("i",Complex(0.,1.)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW), sqrt2("sqrt2",sqrt(2.0));

  Flavour Wbos((kf_code)24);
  Kabbala mw{ "m_{"+Wbos.TexName()+"}", ScalarConstant("mW+")};

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
      if(flav1.Yuk() == 0.0 && flav2.Yuk() == 0.0) continue;

      Kabbala ma{ "m_{"+flav1.TexName()+"}", ScalarConstant("m"+flav1.IDName())};
      Kabbala mb{ "m_{"+flav2.TexName()+"}", ScalarConstant("m"+flav2.IDName())};
      
      /////////////////////////////////////
      // FFS3  -> -PL                    //
      // FFS1  ->  PR                    //
      /////////////////////////////////////
      // // phi +
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav2.Bar());
      m_v.back().AddParticle(flav1);
      m_v.back().AddParticle(Flavour(kf_phiplus));
      m_v.back().Color.push_back
	(i>6?Color_Function(cf::None):
	 Color_Function(cf::D,1,2));
      m_v.back().Color.push_back
	(i>6?Color_Function(cf::None):
	 Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFS3");
      m_v.back().Lorentz.push_back("FFS1");
      m_v.back().cpl.push_back(I*g2*ckm*mb/mw/sqrt2);
      m_v.back().cpl.push_back(-I*g2*ckm*ma/mw/sqrt2);
      m_v.back().order.resize(3);
      m_v.back().order[2] = 1;
      // // phi - 
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav1.Bar());
      m_v.back().AddParticle(flav2);
      m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
      m_v.back().Color.push_back
      	(i>6?Color_Function(cf::None):
      	 Color_Function(cf::D,1,2));
      m_v.back().Color.push_back
      	(i>6?Color_Function(cf::None):
      	 Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFS3");
      m_v.back().Lorentz.push_back("FFS1");
      m_v.back().cpl.push_back(I*g2*ckm*ma/mw/sqrt2);
      m_v.back().cpl.push_back(-I*g2*ckm*mb/mw/sqrt2);
      m_v.back().order.resize(3);
      m_v.back().order[2] = 1;
    }
  }

}


void Standard_ModelGS::GplusPhotonsV()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala two(Kabbala("2",2.0)), I("i",Complex(0.,1.)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW);

  //////////////////////////////////////////////////////////////////////
  // phi+ and phi- vertices with either 1			      //
  // or two photons eqs 75 and 88 https://arxiv.org/pdf/1209.6213.pdf //
  //////////////////////////////////////////////////////////////////////

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(kf_photon);
  m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
  m_v.back().AddParticle(Flavour(kf_phiplus));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VSS2");
  m_v.back().Lorentz.push_back("VSS1");
  m_v.back().cpl.push_back(I*g1);
  m_v.back().cpl.push_back(-I*g1);
  m_v.back().order.resize(3);
  m_v.back().order[2]=1;

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_photon));
  m_v.back().AddParticle(Flavour(kf_photon));
  m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
  m_v.back().AddParticle(Flavour(kf_phiplus));
  m_v.back().cpl.push_back(two*I*g1*g1);
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VVSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]=2;
}

void Standard_ModelGS::GplusZsV()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala one(Kabbala("1",1.0)), two(Kabbala("2",2.0)), I("i",Complex(0.,1.)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW);

  //////////////////////////////////////////////////////////////////////
  // phi+ and phi- vertices with either 1			      //
  // or two Zs eqs 76 and 89 https://arxiv.org/pdf/1209.6213.pdf      //
  //////////////////////////////////////////////////////////////////////
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_Z));
  m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
  m_v.back().AddParticle(Flavour(kf_phiplus));
  m_v.back().cpl.push_back(I*g2/costW*(sqr(costW) - one/two));
  m_v.back().cpl.push_back(-I*g2/costW*(sqr(costW) - one/two));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VSS2");
  m_v.back().Lorentz.push_back("VSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 1;

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_Z));
  m_v.back().AddParticle(Flavour(kf_Z));
  m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
  m_v.back().AddParticle(Flavour(kf_phiplus));
  m_v.back().cpl.push_back(I*g1*g1/two*sqr(costW/sintW-sintW/costW));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VVSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;

  if(Flavour(kf_photon).IsOn()) {
    // // (EQ 95)
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(kf_Z));
    m_v.back().AddParticle(Flavour(kf_phiplus).Bar());
    m_v.back().AddParticle(Flavour(kf_phiplus));
    m_v.back().cpl.push_back( sqr(g1)*I*(costW/sintW - sintW/costW));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 2;
  }
  
}

void Standard_ModelGS::GplusHiggs()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala two(Kabbala("2",2.0)), three(Kabbala("3",3.0));
  Kabbala I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0));
  Kabbala g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED")));
  Kabbala sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW")));
  Kabbala costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW")));
  Kabbala g2(g1/sintW), vev("v_{EW}",ComplexConstant("cvev"));
  Kabbala M{ "M_H", ScalarConstant("mh0") };

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( Flavour(kf_phiplus).Bar() );
  m_v.back().AddParticle( Flavour(kf_phiplus),0);
  m_v.back().AddParticle( Flavour(kf_h0),0) ;
  m_v.back().cpl.push_back(-I*M*M/vev);
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 1;

  // Quadruple phi phi -> h h (EQ 100)
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_phiplus).Bar() );
  m_v.back().AddParticle(Flavour(kf_phiplus) );
  m_v.back().AddParticle(Flavour(kf_h0) );
  m_v.back().AddParticle(Flavour(kf_h0) );
  m_v.back().cpl.push_back( -I*M*M/vev/vev );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;

  // phi h -> W + *
  if(Flavour(kf_Wplus).IsOn()){
    // (EQ 91)
    if(Flavour(kf_Z).IsOn()){
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_h0 ));
      m_v.back().cpl.push_back( -g1*g1/two/costW );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;

      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_h0) ));
      m_v.back().cpl.push_back( g1*g1/two/costW );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;
    }
    // (EQ 93)
    if(Flavour(kf_photon).IsOn()){
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_photon) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
      m_v.back().cpl.push_back( g1*g2/two);
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;

      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_photon) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
      m_v.back().cpl.push_back( -g1*g2/two );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;
    }        
  } // end kf_Wplus
}

void Standard_ModelGS::GplusW()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala one(Kabbala("1",1.0)), two(Kabbala("2",2.0)), I("i",Complex(0.,1.)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW),
    MZ( "M_Z", ScalarConstant("mZ" )),
    vev("v_{EW}",ComplexConstant("cvev"));

  if(Flavour(kf_h0).IsOn()) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
    m_v.back().cpl.push_back( g2/two );
    m_v.back().cpl.push_back( -g2/two );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VSS2");
    m_v.back().Lorentz.push_back("VSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;

    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
    m_v.back().cpl.push_back( g2/two );
    m_v.back().cpl.push_back( -g2/two );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VSS2");
    m_v.back().Lorentz.push_back("VSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;

  } //end kf_h0
  if(Flavour(kf_chi).IsOn()){
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
    m_v.back().cpl.push_back( I*g2/two );
    m_v.back().cpl.push_back( - I*g2/two );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VSS2");
    m_v.back().Lorentz.push_back("VSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;

    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus));
    m_v.back().cpl.push_back( -I*g2/two );
    m_v.back().cpl.push_back( I*g2/two );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VSS2");
    m_v.back().Lorentz.push_back("VSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;
  } //end kf_chi
  if(Flavour(kf_photon).IsOn()){
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_photon) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
    m_v.back().cpl.push_back( g1*g2*vev/two );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;

    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_photon) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus));
    m_v.back().cpl.push_back( g1*g2*vev/two );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;

  } //end kf_photon
  if(Flavour(kf_Z).IsOn()){
    // phi- W+ (EQ 81)
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
    m_v.back().cpl.push_back( -g1*MZ*sintW );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;

    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
    m_v.back().cpl.push_back( g1*MZ*sintW );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;
  } // end kf_Z

  // Quadruple W+ W- (EQ 90)
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus));
  m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
  m_v.back().cpl.push_back( I*g2*g2/two );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VVSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;

  // quadruple phi± chi -> W± A/Z
  if(Flavour(kf_chi).IsOn()){
    // (EQ 92)
    if(Flavour(kf_Z).IsOn()){
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
      m_v.back().cpl.push_back( I*g1*g1/two/costW );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;

      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
      m_v.back().cpl.push_back( I*g1*g1/two/costW );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;

    } // end kf_Z
    // (EQ 94)
    if(Flavour(kf_photon).IsOn()){
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_photon) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
      m_v.back().cpl.push_back( -I*g1*g2/two );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;

      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour(kf_photon) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
      m_v.back().cpl.push_back( -I*g1*g2/two );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS1");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 2;
    }
  }
}

void Standard_ModelGS::GplusGold()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala two(Kabbala("2",2.0)),
    I("i",Complex(0.,1.)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW),
    vev("v_{EW}",ComplexConstant("cvev")),
    MH{ "M_H", ScalarConstant("mh0") };

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
  m_v.back().cpl.push_back( -two*two*I*MH*MH/vev/vev );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;
  if(Flavour(kf_chi)){
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus).Bar() );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_phiplus) );
    m_v.back().cpl.push_back( -two*I*MH*MH/vev/vev );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 2;
  } 
}

void Standard_ModelGS::G0Fermions()
{
  DEBUG_FUNC(METHOD);
  Kabbala two(Kabbala("2",2.0)), three(Kabbala("3",3.0)),
    I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW), vev("v_{EW}",ComplexConstant("cvev")),
    M{ "M_H", ScalarConstant("mh0")};

  std::vector<short int> all_fermions = {1,2,3,4,5,6,11,13,15};

  for(auto i_f: all_fermions){
    Flavour flavi((kf_code)i_f);
    if(flavi.Yuk()){
      Kabbala m{ "m_{"+flavi.TexName()+"}", ScalarConstant("m"+flavi.IDName())};
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)i_f,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)i_f,0) );
      m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
      m_v.back().cpl.push_back( m/vev );
      m_v.back().Color.push_back(i_f>6?Color_Function(cf::None):
				 Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFS2");
      m_v.back().order.resize(3);
      m_v.back().order[2]    = 1;
    }
  }
}

void Standard_ModelGS::G0PhotonsV()
{
  DEBUG_FUNC(METHOD);
}

void Standard_ModelGS::G0ZsV()
{
  DEBUG_FUNC(METHOD);
  Kabbala two(Kabbala("2",2.0)), three(Kabbala("3",3.0)),
    I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW), vev("v_{EW}",ComplexConstant("cvev")),
    M{ "M_H", ScalarConstant("mh0") };

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().cpl.push_back( I*g2*g2/two/costW/costW );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VVSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;
}

void Standard_ModelGS::G0Higgs()
{
  DEBUG_FUNC(METHOD);
  Kabbala two(Kabbala("2",2.0)), three(Kabbala("3",3.0)),
    I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW), vev("v_{EW}",ComplexConstant("cvev")),
    MH{ "M_H", ScalarConstant("mh0") };

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
  m_v.back().cpl.push_back(-I*MH*MH/vev );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 1;

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
  m_v.back().cpl.push_back( -I*MH*MH/vev/vev );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;
  
  if(Flavour(kf_Z).IsOn()){
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle( ATOOLS::Flavour(kf_Z) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
    m_v.back().AddParticle( ATOOLS::Flavour(kf_h0) );
    m_v.back().cpl.push_back( g2/two/costW );
    m_v.back().cpl.push_back( -g2/two/costW );
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VSS2");
    m_v.back().Lorentz.push_back("VSS1");
    m_v.back().order.resize(3);
    m_v.back().order[2]    = 1;
  }
}

void Standard_ModelGS::G0W()
{
  DEBUG_FUNC(METHOD);
  // g1 -> e, g2-> e/sintW
  Kabbala one(Kabbala("1",1.0)), two(Kabbala("2",2.0)), I("i",Complex(0.,1.)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW);

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus).Bar() );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_Wplus) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().cpl.push_back( g2*g2*I/two );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("VVSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;
}

void Standard_ModelGS::G0Gold()
{
  DEBUG_FUNC(METHOD);
  Kabbala two(Kabbala("2",2.0)), three(Kabbala("3",3.0)),
    I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0)),
    g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED"))),
    sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW"))),
    costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW"))),
    g2(g1/sintW), vev("v_{EW}",ComplexConstant("cvev")),
    MH{ "M_H", ScalarConstant("mh0") };
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().AddParticle( ATOOLS::Flavour(kf_chi) );
  m_v.back().cpl.push_back( -three*I*MH*MH/vev/vev );
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSSS1");
  m_v.back().order.resize(3);
  m_v.back().order[2]    = 2;
}

void Standard_ModelGS::ResetVerticesWithEWParameters(const EWParameters& params)
{
  ClearInteractionModel();
  /// Set parameters to their run value before re-initing the vertices
  (*p_complexconstants)[(std::string)( "csin2_thetaW")] = params.m_sw2_r;
  (*p_complexconstants)[(std::string) "ccos2_thetaW"] = params.m_cw2_r;
  (*p_complexconstants)[(std::string) "cvev"]         = params.m_cvev_r;
  (*p_constants)[(std::string) ("mZ")]                  = params.m_mz_r;
  (*p_constants)[(std::string) "mW+"]                 = params.m_mw_r;
  (*p_constants)[(std::string) "mh0"]                 = params.m_mh0_r;
  (*p_constants)[(std::string) "mt"]                  = params.m_mt_r;
  (*p_constants)[(std::string) "alpha_QED"]           = params.m_aew_r;
  InitializeInteractionModel();
}

void Standard_ModelGS::ClearInteractionModel()
{
  m_v.clear();
  m_ov.clear();
  m_fls.clear();
  m_vmap.clear();
  m_vtable.clear();
}

size_t Standard_ModelGS::IndexOfOrderKey(const std::string& key) const
{
  PRINT_VAR(key);
  if(key == "SMGold")
    return 2;
  else return Model_Base::IndexOfOrderKey(key);
}
