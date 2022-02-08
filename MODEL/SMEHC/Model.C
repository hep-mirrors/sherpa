#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace MODEL {

  class SMEHC: public Model_Base {
  private:

    int  m_ckmorder;

    void RegisterDefaults() const override;

    void FixEWParameters();
    void FixCKM();
    void FixEFT();

    void ParticleInit() override;

    void InitQEDVertices();
    void InitQCDVertices();
    void InitEWVertices();
    void InitEFTVertices();

  public :

    SMEHC();
    bool ModelInit() override;
    void InitVertices() override;
    size_t IndexOfOrderKey(const std::string&) const override;

  };

}

#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "MODEL/SMEHC/Effective_Higgs_Coupling.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(SMEHC,"SMEHC",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,SMEHC>::
operator()(const Model_Arguments &args) const
{
  return new SMEHC();
}

void Getter<Model_Base,Model_Arguments,SMEHC>::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"The Standard Model\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"# possible parameters in yaml configuration [usage: \"keyword: value\"]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (EW input scheme, see documentation)\n"
     <<setw(width+7)<<" "<<"- EW_REN_SCHEME (EW renormalisation scheme, see documentation)\n"
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

SMEHC::SMEHC() :
  Model_Base(1)
{
  m_name="SMEHC";
  ParticleInit();
  Model_Base::RegisterDefaults();
  AddStandardContainers();
  CustomContainerInit();
}

void SMEHC::RegisterDefaults() const
{
  Settings& s = Settings::GetMainSettings();
  s["FINITE_TOP_MASS"].SetDefault(false);
  s["FINITE_W_MASS"].SetDefault(false);
  s["DEACTIVATE_GGH"].SetDefault(false);
  s["DEACTIVATE_PPH"].SetDefault(false);
  s["ALPHAS_GGH"].SetDefault(ScalarConstant("alpha_S"));
  s["1/ALPHAQED_PPH"].SetDefault(1.0/ScalarConstant("alpha_QED"));
}

void SMEHC::ParticleInit()
{
  s_kftable[kf_none] = new ATOOLS::Particle_Info(kf_none,-1,0,0,0,0,-1,0,1,0,"no_particle","no_particle","no_particle", "no_particle", 1,1);
  //add SM particles
  //kf_code,mass,width,charge,strong,spin,majorana,take,stable,massive,idname,antiname,texname,antitexname
  s_kftable[kf_d]      = new Particle_Info(kf_d,0.01,.0,-1,3,1,0,1,1,0,"d","db", "d", "\\bar{d}");
  s_kftable[kf_u]      = new Particle_Info(kf_u,0.005,.0,2,3,1,0,1,1,0,"u","ub", "u", "\\bar{u}");
  s_kftable[kf_s]      = new Particle_Info(kf_s,0.2,.0,-1,3,1,0,1,1,0,"s","sb", "s", "\\bar{s}");
  s_kftable[kf_c]      = new Particle_Info(kf_c,1.42,.0,2,3,1,0,1,1,0,"c","cb", "c", "\\bar{c}");
  s_kftable[kf_b]      = new Particle_Info(kf_b,4.92,.0,-1,3,1,0,1,1,0,"b","bb", "b", "\\bar{b}");
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
  s_kftable[kf_shgluon] = new Particle_Info(kf_shgluon,0.0,0.0,0,8,2,-1,1,1,0,"GS","GS","G_{S}","G_{S}",1);
  s_kftable[kf_gluon_qgc] = new Particle_Info(kf_gluon_qgc,0.0,0.0,0,8,4,-1,1,1,0,"G4","G4","G_{4}","G_{4}",1);
  ReadParticleData();
}

bool SMEHC::ModelInit()
{
  Settings& s = Settings::GetMainSettings();
  FixEWParameters();  
  FixCKM();
  SetAlphaQCD(*p_isrhandlermap, s["ALPHAS(MZ)"].Get<double>());
  RegisterDefaults();
  FixEFT();
  SetRunningFermionMasses();
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

void SMEHC::FixEWParameters()
{
  Settings& s = Settings::GetMainSettings();
  Complex csin2thetaW, ccos2thetaW, cvev, I(0.,1.);
  string yukscheme = s["YUKAWA_MASSES"].Get<string>();
  p_numbers->insert(make_pair(string("YukawaScheme"), yukscheme=="Running"));
  string widthscheme = s["WIDTH_SCHEME"].Get<string>();
  p_numbers->insert(make_pair(string("WidthScheme"), widthscheme=="CMS"));
  ew_scheme::code ewscheme = s["EW_SCHEME"].Get<ew_scheme::code>();
  ew_scheme::code ewrenscheme = s["EW_REN_SCHEME"].Get<ew_scheme::code>();
  double MW=Flavour(kf_Wplus).Mass(), GW=Flavour(kf_Wplus).Width();
  double MZ=Flavour(kf_Z).Mass(), GZ=Flavour(kf_Z).Width();
  double MH=Flavour(kf_h0).Mass(), GH=Flavour(kf_h0).Width();
  std::string ewschemename(""),ewrenschemename("");
  PRINT_VAR(ewscheme);
  switch (ewscheme) {
  case ew_scheme::UserDefined:
    // all SM parameters given explicitly
    ewschemename="user-defined, input: all parameters";
    SetAlphaQEDByScale(s["ALPHAQED_DEFAULT_SCALE"].Get<double>());
    csin2thetaW = s["SIN2THETAW"].Get<double>();
    ccos2thetaW=1.-csin2thetaW;
    cvev = s["VEV"].Get<double>();
    break;
  case ew_scheme::alpha0: {
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
  case ew_scheme::alphamZ: {
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
  case ew_scheme::Gmu: {
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
  case ew_scheme::alphamZsW: {
    // alpha(mZ)-mZ-sin(theta) scheme
    ewschemename="alpha(mZ)-mZ-sin(theta_W) scheme, input: 1/\\alphaQED(m_Z), sin^2(theta_W), m_Z, m_h, widths";
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
  case ew_scheme::alphamWsW: {
    // alpha(mW)-mW-sin(theta) scheme
    ewschemename="alpha(mW)-mW-sin(theta_W) scheme, input: 1/\\alphaQED(m_W), sin^2(theta_W), m_W, m_h, widths";
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
  case ew_scheme::GmumZsW: {
    // Gmu-mZ-sin(theta) scheme
    ewschemename="Gmu-mZ-sin(theta_W) scheme, input: GF, sin^2(theta_W), m_Z, m_h, widths";
    double GF = s["GF"].Get<double>();
    csin2thetaW = s["SIN2THETAW"].Get<double>();
    ccos2thetaW=1.-csin2thetaW;
    MW=MZ*sqrt(ccos2thetaW.real());
    Flavour(kf_Wplus).SetMass(MW);
    cvev=1./(pow(2.,0.25)*sqrt(GF));
    if (widthscheme=="CMS") {
      Complex muW2(0.,0.), muZ2(MZ*(MZ-I*GZ));
      muW2=muZ2*ccos2thetaW;
      MW=sqrt(muW2.real());
      GW=-muW2.imag()/MW;
      Flavour(kf_Wplus).SetMass(MW);
      Flavour(kf_Wplus).SetWidth(GW);
      cvev=1./(pow(2.,0.25)*sqrt(GF));
      const size_t aqedconv{ s["GMU_CMS_AQED_CONVENTION"].Get<size_t>() };
      switch (aqedconv) {
      case 0:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::abs(muZ2*csin2thetaW*(1.-csin2thetaW)));
        break;
      case 1:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muZ2*csin2thetaW*(1.-csin2thetaW)));
        break;
      case 2:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muZ2*(1.-csin2thetaW))*std::real(csin2thetaW));
        break;
      case 3 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MZ)*(1.-csin2thetaW.real())*std::abs(csin2thetaW));
        break;
      case 4 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MZ)*(1.-csin2thetaW.real())*csin2thetaW.real());
        break;
      case 5:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muZ2)*std::real(1.-csin2thetaW)*std::real(csin2thetaW));
        break;
      case 6 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MZ)*std::abs((1.-csin2thetaW)*csin2thetaW));
        break;
      default:
        THROW(not_implemented,"\\alpha_QED convention not implemented.");
      }
    } else if (widthscheme=="Fixed") {
      if (csin2thetaW.imag()!=0.0) THROW(fatal_error,"sin^2(\\theta_w) not real.");
      SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MZ)*csin2thetaW.real()*(1.-csin2thetaW.real()));
    }
    break;
  }
  case ew_scheme::GmumWsW: {
    // Gmu-mW-sin(theta) scheme
    ewschemename="Gmu-mW-sin(theta_W) scheme, input: GF, sin^2(theta_W), m_W, m_h, widths";
    double GF = s["GF"].Get<double>();
    csin2thetaW = s["SIN2THETAW"].Get<double>();
    ccos2thetaW=1.-csin2thetaW;
    MZ=MW/sqrt(ccos2thetaW.real());
    Flavour(kf_Z).SetMass(MZ);
    cvev=1./(pow(2.,0.25)*sqrt(GF));
    if (widthscheme=="CMS") {
      Complex muW2(MW*(MW-I*GW)), muZ2(0.,0.);
      muZ2=muW2/ccos2thetaW;
      MZ=sqrt(muZ2.real());
      GZ=-muZ2.imag()/MZ;
      Flavour(kf_Z).SetMass(MZ);
      Flavour(kf_Z).SetWidth(GZ);
      cvev=1./(pow(2.,0.25)*sqrt(GF));
      const size_t aqedconv{ s["GMU_CMS_AQED_CONVENTION"].Get<size_t>() };
      switch (aqedconv) {
      case 0:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::abs(muW2)*csin2thetaW.real());
        break;
      case 1:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muW2)*csin2thetaW.real());
        break;
      case 2:
        SetAlphaQED(sqrt(2.)*GF/M_PI*std::real(muW2)*csin2thetaW.real());
        break;
      case 3 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MW)*csin2thetaW.real());
        break;
      case 4 :
        SetAlphaQED(sqrt(2.)*GF/M_PI*sqr(MW)*csin2thetaW.real());
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
  case ew_scheme::FeynRules: {
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
    msg_Info()<<"Unknown EW_REN_SCHEME="<<ewrenscheme<<", resetting to Gmu."
              <<std::endl;
    ewrenscheme=ew_scheme::Gmu;
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

void SMEHC::FixCKM()
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

void SMEHC::FixEFT()
{
  Settings& s = Settings::GetMainSettings();
  Complex eh(2./3.,0.);
  double ph(-2.*47./18.); double pfac(1.);
  if (s["FINITE_TOP_MASS"].Get<bool>()) {
    double hm=Flavour(kf_h0).Mass();
    double tm=Flavour(kf_t).Yuk();
    Effective_Higgs_Coupling ehc(hm);
    eh = ehc.GetFermionContribution(tm);
    double taut(sqr(hm/2./tm));
    pfac += -56./705.*taut - 32./987.*sqr(taut);
  }
  if (s["FINITE_W_MASS"].Get<bool>()) {
    double hm=Flavour(kf_h0).Mass();
    double Wm=Flavour(kf_Wplus).Mass();
    double tauW(sqr(hm/2./Wm));
    pfac += 66./235.*tauW + 228./1645.*sqr(tauW) + 696./8225.*tauW*sqr(tauW)
            + 5248./90475.*sqr(sqr(tauW)) + 1280./29939.*tauW*sqr(sqr(tauW))
            + 54528./1646645.*sqr(tauW*sqr(tauW));
  }
  if (s["DEACTIVATE_GGH"].Get<int>())
    eh=Complex(0.,0.);
  if (s["DEACTIVATE_PPH"].Get<int>())
    ph=0.;
  p_constants->insert(std::make_pair(std::string("h0_pp_fac"),ph*pfac));
  p_constants->insert(std::make_pair(std::string("h0_gg_fac"),real(eh)));
}

void SMEHC::InitVertices()
{
  InitQEDVertices();
  InitQCDVertices();
  InitEWVertices();
  InitEFTVertices();
}

void SMEHC::InitQEDVertices()
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

void SMEHC::InitQCDVertices()
{
  if (!Flavour(kf_gluon).IsOn()) return;
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
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<2;++i) m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_gluon_qgc));
  m_v.back().Color.push_back(Color_Function(cf::F,1,2,3));
  m_v.back().Lorentz.push_back("VVP");
  m_v.back().cpl.push_back(cpl1);
  m_v.back().order[0]=1;
  m_v.back().dec=1;
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
  m_v.back().dec=-1;
}

void SMEHC::InitEWVertices()
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

void SMEHC::InitEFTVertices()
{
  if (!Flavour(kf_h0).IsOn()) return;
  Settings& s = Settings::GetMainSettings();
  const double as{ s["ALPHAS_GGH"].Get<double>() };
  const double aqed{ 1.0/s["1/ALPHAQED_PPH"].Get<double>() };
  Kabbala ghgg("g_hgg",ScalarConstant("h0_gg_fac")*
	       as/(2.0*M_PI*ComplexConstant("cvev")));
  Kabbala ghpp("g_hpp",ScalarConstant("h0_pp_fac")*
	       aqed/(2.0*M_PI*ComplexConstant("cvev")));
  Kabbala g3("g_3",sqrt(4.*M_PI*ScalarConstant("alpha_S")));
  Kabbala I("i",Complex(0.,1.));
  msg_Info()<<METHOD<<"() {\n  ggh coupling is "
	    <<ghgg.Value()<<" [ \\alpha_s = "<<as
	    <<" ]\n  yyh coupling is "<<ghpp.Value()
            <<" [ 1/\\alpha = "<<(1.0/aqed)<<" ]\n}\n";
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_h0));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().Color.push_back(Color_Function(cf::G,2,3));
  m_v.back().Lorentz.push_back("HVV");
  m_v.back().cpl.push_back(I*ghgg);
  m_v.back().order[0]=2;
  m_v.back().order[1]=1;
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_h0));
  m_v.back().AddParticle(Flavour(kf_photon));
  m_v.back().AddParticle(Flavour(kf_photon));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("HVV");
  m_v.back().cpl.push_back(I*ghpp);
  m_v.back().order[1]=3;
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_h0));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().Color.push_back(Color_Function(cf::F,2,3,4));
  m_v.back().Lorentz.push_back("HVVV");
  m_v.back().cpl.push_back(ghgg*g3);
  m_v.back().order[0]=3;
  m_v.back().order[1]=1;
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_h0));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_gluon_qgc));
  m_v.back().Color.push_back(Color_Function(cf::F,2,3,4));
  m_v.back().Lorentz.push_back("HVVP");
  m_v.back().cpl.push_back(ghgg*g3);
  m_v.back().order[0]=3;
  m_v.back().order[1]=1;
  m_v.back().dec=3;
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<4;++i) m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_h0));
  for (size_t i(0);i<3;++i) m_v.back().cpl.push_back(-I*ghgg*g3*g3);
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,2,new Color_Function(cf::F,3,4,-1)));
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,3,new Color_Function(cf::F,2,4,-1)));
  m_v.back().Color.push_back
    (Color_Function(cf::F,-1,1,4,new Color_Function(cf::F,2,3,-1)));
  m_v.back().Lorentz.push_back("HVVVVA");
  m_v.back().Lorentz.push_back("HVVVVB");
  m_v.back().Lorentz.push_back("HVVVVC");
  m_v.back().order[0]=4;
  m_v.back().order[1]=1;
  m_v.back().dec=3;
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_h0));
  m_v.back().AddParticle(Flavour(kf_shgluon));
  m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().Color.push_back(Color_Function(cf::G,2,3));
  m_v.back().Lorentz.push_back("C4GS");
  m_v.back().cpl.push_back(-I);
  m_v.back().dec=-1;
  Kabbala cpl2=-I*g3*g3*ghgg;
  m_v.push_back(Single_Vertex());
  for (size_t i(0);i<3;++i) m_v.back().AddParticle(Flavour(kf_gluon));
  m_v.back().AddParticle(Flavour(kf_shgluon));
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
  m_v.back().order[1]=1;
  m_v.back().dec=-1;
}

size_t SMEHC::IndexOfOrderKey(const std::string& key) const
{
  if (key == "SMEHC")
    return 2;
  else
    return Model_Base::IndexOfOrderKey(key);
}
