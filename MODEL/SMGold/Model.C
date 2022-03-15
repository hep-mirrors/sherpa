#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
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

#include "MODEL/Main/Single_Vertex.H"
#include <iomanip>

using namespace ATOOLS;
using namespace std;

namespace MODEL{
  using complex = std::complex<double>;
  
  class Standard_ModelGS : public Model_Base
  {
  private:
    ComplexConstantsMap m_local_complex_constants;
    int  m_ckmorder, m_dec_g4;
  public:
  Standard_ModelGS() : Model_Base(true)
      {
	m_name = string("SMGold");
	ParticleInit();
	RegisterDefaults();
	// Massive and Stable flags
	// are set consistently with
	// UFO above, ReadParticleData
	// allows to overwrite these flags.
	// Need this before AddStandardContainers.
	AddStandardContainers();
	CustomContainerInit();
      }
  protected:
    void ParticleInit()
    {
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
      s_kftable[kf_t]         = new Particle_Info(kf_t,172.5,1.32,2,3,1,0,1,0,1,"t","tb", "t", "\\bar{t}");
      s_kftable[kf_e]         = new Particle_Info(kf_e,0.000511,.0,-3,0,1,0,1,1,0,"e-","e+", "e^{-}", "e^{+}");
      s_kftable[kf_nue]       = new Particle_Info(kf_nue,.0,.0,0,0,1,0,1,1,0,"ve","veb", "\\nu_{e}", "\\bar{\\nu}_{e}");
      s_kftable[kf_mu]        = new Particle_Info(kf_mu,.105,.0,-3,0,1,0,1,1,0,"mu-","mu+", "\\mu^{-}", "\\mu^{+}");
      s_kftable[kf_numu]      = new Particle_Info(kf_numu,.0,.0,0,0,1,0,1,1,0,"vmu","vmub", "\\nu_{\\mu}", "\\bar{\\nu}_{\\mu}");
      s_kftable[kf_tau]       = new Particle_Info(kf_tau,1.777,2.26735e-12,-3,0,1,0,1,0,0,"tau-","tau+", "\\tau^{-}", "\\tau^{+}");
      s_kftable[kf_nutau]     = new Particle_Info(kf_nutau,.0,.0,0,0,1,0,1,1,0,"vtau","vtaub", "\\nu_{\\tau}", "\\bar{\\nu}_{\\tau}");
      s_kftable[kf_gluon]     = new Particle_Info(kf_gluon,.0,.0,0,8,2,-1,1,1,0,"G","G", "G", "G");
      s_kftable[kf_photon]    = new Particle_Info(kf_photon,.0,.0,0,0,2,-1,1,1,0,"P","P","\\gamma","\\gamma");
      s_kftable[kf_Z]         = new Particle_Info(kf_Z,91.1876,2.4952,0,0,2,-1,1,0,1,"Z","Z","Z","Z");
      s_kftable[kf_Wplus]     = new Particle_Info(kf_Wplus,80.379,2.085,3,0,2,0,1,0,1,"W+","W-","W^{+}","W^{-}");
      s_kftable[kf_gluon_qgc] = new Particle_Info(kf_gluon_qgc,0.0,0.0,0,8,4,-1,1,1,0,"G4","G4","G_{4}","G_{4}",1);
  
      s_kftable[kf_phiplus]   = new Particle_Info(kf_phiplus,80.385,2.085,3,0,0,0,1,0,1,"phi+","phi-","\\phi^{+}","\\phi^{-}");
      s_kftable[kf_chi]       = new Particle_Info(kf_chi,91.1876,2.4952,0,0,0,-1,1,0,1,"chi","chi","\\chi","\\chi");
      ReadParticleData();
    }
    bool ModelInit()
    {
      FixEWParameters();  
      FixCKM();
      Settings& s = Settings::GetMainSettings();
      SetAlphaQCD(*p_isrhandlermap, s["ALPHAS(MZ)"].Get<double>());
      SetRunningFermionMasses();
      SetRunningBosonMasses();
      ATOOLS::OutputParticles(msg->Info());
      ATOOLS::OutputContainers(msg->Info());
      OutputCKM();
      //      ParamInit();
      for (MODEL::ScalarNumbersMap::iterator it=p_numbers->begin();
	   it!=p_numbers->end();++it) DEBUG_INFO(it->first+" = "<<it->second);
      for (MODEL::ScalarConstantsMap::iterator it=p_constants->begin();
	   it!=p_constants->end();++it) DEBUG_INFO(it->first+" = "<<it->second);
      for (MODEL::ComplexConstantsMap::iterator it=p_complexconstants->begin();
	   it!=p_complexconstants->end();++it) DEBUG_INFO(it->first+" = "<<it->second);
      return true;
    }

    void FixEWParameters()
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
    void FixCKM()
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
    void ParamInit()
    {
      DEBUG_FUNC(this);
      m_local_complex_constants.clear();
      msg_Debugging() << setprecision(20);
      double ymb = Flavour(kf_b).Yuk();
      double ymt = ScalarConstant("mt");
      double MH = ScalarConstant("mh0");
      Complex CKM1x1 = ComplexConstant("CKM_1_1");
      DEBUG_VAR(CKM1x1);
      Complex CKM1x2 = ComplexConstant("CKM_1_2");
      DEBUG_VAR(CKM1x2);
      Complex CKM2x1 = ComplexConstant("CKM_2_1");
      DEBUG_VAR(CKM2x1);
      Complex CKM2x2 = ComplexConstant("CKM_2_2");
      DEBUG_VAR(CKM2x2);
      double aEW = ScalarConstant("alpha_QED");
      DEBUG_VAR(aEW);
      double G = (((2.0*sqrt(ScalarConstant("alpha_S")))*sqrt(M_PI)));
      DEBUG_VAR(G);
      double ee = (((2.0*sqrt(aEW))*sqrt(M_PI)));
      DEBUG_VAR(ee);
      double sw2 = ComplexConstant("csin2_thetaW").real();
      DEBUG_VAR(sw2);
      double cw2 = ComplexConstant("ccos2_thetaW").real();;
      double cw  = sqrt(cw2);
      DEBUG_VAR(cw);
      double sw = (sqrt(sw2));
      DEBUG_VAR(sw);
      double g1 = ((ee/cw));
      DEBUG_VAR(g1);
      double gw = ((ee/sw));
      DEBUG_VAR(gw);
      double vev = ComplexConstant("cvev").real();
      DEBUG_VAR(vev);
      double lam = ((pow(MH,2.0)/(2.0*pow(vev,2.0))));
      DEBUG_VAR(lam);
      double yb = (((ymb*sqrt(2.0))/vev));
      DEBUG_VAR(yb);
      double yt = (((ymt*sqrt(2.0))/vev));
      DEBUG_VAR(yt);
      double muH = (sqrt((lam*pow(vev,2.0))));
      DEBUG_VAR(muH);
      Complex I1a33 = yb;
      DEBUG_VAR(I1a33);
      Complex I2a33 = yt;
      DEBUG_VAR(I2a33);
      Complex I3a33 = yt;
      DEBUG_VAR(I3a33);
      Complex I4a33 = yb;
      DEBUG_VAR(I4a33);
      
      m_local_complex_constants.insert(make_pair(string("GC_22"),((-6.0*complex(0.0,1.0))*lam)));
      DEBUG_VAR((m_local_complex_constants)["GC_22"]);
      m_local_complex_constants.insert(make_pair(string("GC_20"),((-2.0*complex(0.0,1.0))*lam)));
      DEBUG_VAR((m_local_complex_constants)["GC_20"]);
      m_local_complex_constants.insert(make_pair(string("GC_21"),((-4.0*complex(0.0,1.0))*lam)));
      DEBUG_VAR((m_local_complex_constants)["GC_21"]);
      m_local_complex_constants.insert(make_pair(string("GC_20"),((-2.0*complex(0.0,1.0))*lam)));
      DEBUG_VAR((m_local_complex_constants)["GC_20"]);
      m_local_complex_constants.insert(make_pair(string("GC_20"),((-2.0*complex(0.0,1.0))*lam)));
      DEBUG_VAR((m_local_complex_constants)["GC_20"]);
      m_local_complex_constants.insert(make_pair(string("GC_22"),((-6.0*complex(0.0,1.0))*lam)));
      DEBUG_VAR((m_local_complex_constants)["GC_22"]);
      m_local_complex_constants.insert(make_pair(string("GC_59"),(((-2.0*complex(0.0,1.0))*lam)*vev)));
      DEBUG_VAR((m_local_complex_constants)["GC_59"]);
      m_local_complex_constants.insert(make_pair(string("GC_59"),(((-2.0*complex(0.0,1.0))*lam)*vev)));
      DEBUG_VAR((m_local_complex_constants)["GC_59"]);
      m_local_complex_constants.insert(make_pair(string("GC_60"),(((-6.0*complex(0.0,1.0))*lam)*vev)));
      DEBUG_VAR((m_local_complex_constants)["GC_60"]);
      m_local_complex_constants.insert(make_pair(string("GC_7"),((2.0*pow(ee,2.0))*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_7"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_66"),((pow(ee,2.0)*vev)/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_66"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_61"),((-(pow(ee,2.0)*vev))/(4.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_61"]);
      m_local_complex_constants.insert(make_pair(string("GC_62"),((-((pow(ee,2.0)*complex(0.0,1.0))*vev))/(4.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_62"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_69"),(((-(pow(ee,2.0)*vev))/(4.0*cw))+(((cw*pow(ee,2.0))*vev)/(4.0*pow(sw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_69"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_65"),((-(pow(ee,2.0)*vev))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_65"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_64"),((pow(ee,2.0)*vev)/(4.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_64"]);
      m_local_complex_constants.insert(make_pair(string("GC_62"),((-((pow(ee,2.0)*complex(0.0,1.0))*vev))/(4.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_62"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_68"),(((pow(ee,2.0)*vev)/(4.0*cw))-(((cw*pow(ee,2.0))*vev)/(4.0*pow(sw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_68"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_70"),(((pow(ee,2.0)*vev)/(4.0*cw))+(((cw*pow(ee,2.0))*vev)/(4.0*pow(sw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_70"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_67"),(((-(pow(ee,2.0)*vev))/(4.0*cw))-(((cw*pow(ee,2.0))*vev)/(4.0*pow(sw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_67"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_71"),((((-((pow(ee,2.0)*complex(0.0,1.0))*vev))/2.0)-((((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))*vev)/(4.0*pow(sw,2.0))))-((((pow(ee,2.0)*complex(0.0,1.0))*pow(sw,2.0))*vev)/(4.0*pow(cw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_71"]);
      m_local_complex_constants.insert(make_pair(string("GC_11"),(-G)));
      DEBUG_VAR((m_local_complex_constants)["GC_11"]);
      m_local_complex_constants.insert(make_pair(string("GC_11"),(-G)));
      DEBUG_VAR((m_local_complex_constants)["GC_11"]);
      m_local_complex_constants.insert(make_pair(string("GC_13"),G));
      DEBUG_VAR((m_local_complex_constants)["GC_13"]);
      m_local_complex_constants.insert(make_pair(string("GC_11"),(-G)));
      DEBUG_VAR((m_local_complex_constants)["GC_11"]);
      m_local_complex_constants.insert(make_pair(string("GC_13"),G));
      DEBUG_VAR((m_local_complex_constants)["GC_13"]);
      m_local_complex_constants.insert(make_pair(string("GC_11"),(-G)));
      DEBUG_VAR((m_local_complex_constants)["GC_11"]);
      m_local_complex_constants.insert(make_pair(string("GC_11"),(-G)));
      DEBUG_VAR((m_local_complex_constants)["GC_11"]);
      m_local_complex_constants.insert(make_pair(string("GC_13"),G));
      DEBUG_VAR((m_local_complex_constants)["GC_13"]);
      m_local_complex_constants.insert(make_pair(string("GC_14"),(-(complex(0.0,1.0)*pow(G,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_14"]);
      m_local_complex_constants.insert(make_pair(string("GC_14"),(-(complex(0.0,1.0)*pow(G,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_14"]);
      m_local_complex_constants.insert(make_pair(string("GC_15"),(complex(0.0,1.0)*pow(G,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_15"]);
      m_local_complex_constants.insert(make_pair(string("GC_15"),(complex(0.0,1.0)*pow(G,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_15"]);
      m_local_complex_constants.insert(make_pair(string("GC_14"),(-(complex(0.0,1.0)*pow(G,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_14"]);
      m_local_complex_constants.insert(make_pair(string("GC_15"),(complex(0.0,1.0)*pow(G,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_15"]);
      m_local_complex_constants.insert(make_pair(string("GC_17"),(-I2a33)));
      DEBUG_VAR((m_local_complex_constants)["GC_17"]);
      m_local_complex_constants.insert(make_pair(string("GC_16"),I1a33));
      DEBUG_VAR((m_local_complex_constants)["GC_16"]);
      m_local_complex_constants.insert(make_pair(string("GC_75"),(yb/sqrt(2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_75"]);
      m_local_complex_constants.insert(make_pair(string("GC_73"),(-(yb/sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_73"]);
      m_local_complex_constants.insert(make_pair(string("GC_74"),(-((complex(0.0,1.0)*yb)/sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_74"]);
      m_local_complex_constants.insert(make_pair(string("GC_74"),(-((complex(0.0,1.0)*yb)/sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_74"]);
      m_local_complex_constants.insert(make_pair(string("GC_19"),(-I4a33)));
      DEBUG_VAR((m_local_complex_constants)["GC_19"]);
      m_local_complex_constants.insert(make_pair(string("GC_18"),I3a33));
      DEBUG_VAR((m_local_complex_constants)["GC_18"]);
      m_local_complex_constants.insert(make_pair(string("GC_76"),(-(yt/sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_76"]);
      m_local_complex_constants.insert(make_pair(string("GC_78"),(yt/sqrt(2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_78"]);
      m_local_complex_constants.insert(make_pair(string("GC_77"),(-((complex(0.0,1.0)*yt)/sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_77"]);
      m_local_complex_constants.insert(make_pair(string("GC_77"),(-((complex(0.0,1.0)*yt)/sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_77"]);
      m_local_complex_constants.insert(make_pair(string("GC_40"),((-(pow(ee,2.0)*complex(0.0,1.0)))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_40"]);
      m_local_complex_constants.insert(make_pair(string("GC_39"),((-pow(ee,2.0))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_39"]);
      m_local_complex_constants.insert(make_pair(string("GC_65"),((-(pow(ee,2.0)*vev))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_65"]);
      m_local_complex_constants.insert(make_pair(string("GC_29"),((-(ee*complex(0.0,1.0)))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_29"]);
      m_local_complex_constants.insert(make_pair(string("GC_30"),((ee*complex(0.0,1.0))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_30"]);
      m_local_complex_constants.insert(make_pair(string("GC_31"),(ee/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_31"]);
      m_local_complex_constants.insert(make_pair(string("GC_28"),((-ee)/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_28"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_4"),(ee*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_4"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_40"),((-(pow(ee,2.0)*complex(0.0,1.0)))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_40"]);
      m_local_complex_constants.insert(make_pair(string("GC_41"),(pow(ee,2.0)/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_41"]);
      m_local_complex_constants.insert(make_pair(string("GC_66"),((pow(ee,2.0)*vev)/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_66"]);
      m_local_complex_constants.insert(make_pair(string("GC_30"),((ee*complex(0.0,1.0))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_30"]);
      m_local_complex_constants.insert(make_pair(string("GC_29"),((-(ee*complex(0.0,1.0)))/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_29"]);
      m_local_complex_constants.insert(make_pair(string("GC_31"),(ee/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_31"]);
      m_local_complex_constants.insert(make_pair(string("GC_28"),((-ee)/(2.0*sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_28"]);
      m_local_complex_constants.insert(make_pair(string("GC_23"),((pow(ee,2.0)*complex(0.0,1.0))/(2.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_23"]);
      m_local_complex_constants.insert(make_pair(string("GC_23"),((pow(ee,2.0)*complex(0.0,1.0))/(2.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_23"]);
      m_local_complex_constants.insert(make_pair(string("GC_23"),((pow(ee,2.0)*complex(0.0,1.0))/(2.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_23"]);
      m_local_complex_constants.insert(make_pair(string("GC_63"),(((pow(ee,2.0)*complex(0.0,1.0))*vev)/(2.0*pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_63"]);
      m_local_complex_constants.insert(make_pair(string("GC_5"),(pow(ee,2.0)*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_5"]);
      m_local_complex_constants.insert(make_pair(string("GC_5"),(pow(ee,2.0)*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_5"]);
      m_local_complex_constants.insert(make_pair(string("GC_6"),((-2.0*pow(ee,2.0))*complex(0.0,1.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_6"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_38"),(((cw*ee)*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_38"]);
      m_local_complex_constants.insert(make_pair(string("GC_37"),(-(((cw*ee)*complex(0.0,1.0))/sw))));
      DEBUG_VAR((m_local_complex_constants)["GC_37"]);
      m_local_complex_constants.insert(make_pair(string("GC_24"),(-((pow(ee,2.0)*complex(0.0,1.0))/pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_24"]);
      m_local_complex_constants.insert(make_pair(string("GC_24"),(-((pow(ee,2.0)*complex(0.0,1.0))/pow(sw,2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_24"]);
      m_local_complex_constants.insert(make_pair(string("GC_25"),(((2.0*pow(ee,2.0))*complex(0.0,1.0))/pow(sw,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_25"]);
      m_local_complex_constants.insert(make_pair(string("GC_54"),((((cw*pow(ee,2.0))*complex(0.0,1.0))/sw)-(((pow(ee,2.0)*complex(0.0,1.0))*sw)/cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_54"]);
      m_local_complex_constants.insert(make_pair(string("GC_53"),(((cw*ee)/(2.0*sw))+((ee*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_53"]);
      m_local_complex_constants.insert(make_pair(string("GC_47"),(((-(cw*ee))/(2.0*sw))-((ee*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_47"]);
      m_local_complex_constants.insert(make_pair(string("GC_50"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_50"]);
      m_local_complex_constants.insert(make_pair(string("GC_51"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_51"]);
      m_local_complex_constants.insert(make_pair(string("GC_9"),((pow(ee,2.0)*complex(0.0,1.0))/(2.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_9"]);
      m_local_complex_constants.insert(make_pair(string("GC_10"),(pow(ee,2.0)/(2.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_10"]);
      m_local_complex_constants.insert(make_pair(string("GC_58"),((pow(ee,2.0)*vev)/(2.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_58"]);
      m_local_complex_constants.insert(make_pair(string("GC_9"),((pow(ee,2.0)*complex(0.0,1.0))/(2.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_9"]);
      m_local_complex_constants.insert(make_pair(string("GC_8"),((-pow(ee,2.0))/(2.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_8"]);
      m_local_complex_constants.insert(make_pair(string("GC_57"),((-(pow(ee,2.0)*vev))/(2.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_57"]);
      m_local_complex_constants.insert(make_pair(string("GC_42"),(((cw*pow(ee,2.0))*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_42"]);
      m_local_complex_constants.insert(make_pair(string("GC_43"),((((-2.0*cw)*pow(ee,2.0))*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_43"]);
      m_local_complex_constants.insert(make_pair(string("GC_42"),(((cw*pow(ee,2.0))*complex(0.0,1.0))/sw)));
      DEBUG_VAR((m_local_complex_constants)["GC_42"]);
      m_local_complex_constants.insert(make_pair(string("GC_56"),(((pow(ee,2.0)*complex(0.0,1.0))+(((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))/(2.0*pow(sw,2.0))))+(((pow(ee,2.0)*complex(0.0,1.0))*pow(sw,2.0))/(2.0*pow(cw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_56"]);
      m_local_complex_constants.insert(make_pair(string("GC_55"),(((-(pow(ee,2.0)*complex(0.0,1.0)))+(((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))/(2.0*pow(sw,2.0))))+(((pow(ee,2.0)*complex(0.0,1.0))*pow(sw,2.0))/(2.0*pow(cw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_55"]);
      m_local_complex_constants.insert(make_pair(string("GC_56"),(((pow(ee,2.0)*complex(0.0,1.0))+(((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))/(2.0*pow(sw,2.0))))+(((pow(ee,2.0)*complex(0.0,1.0))*pow(sw,2.0))/(2.0*pow(cw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_56"]);
      m_local_complex_constants.insert(make_pair(string("GC_72"),((((pow(ee,2.0)*complex(0.0,1.0))*vev)+((((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))*vev)/(2.0*pow(sw,2.0))))+((((pow(ee,2.0)*complex(0.0,1.0))*pow(sw,2.0))*vev)/(2.0*pow(cw,2.0))))));
      DEBUG_VAR((m_local_complex_constants)["GC_72"]);
      m_local_complex_constants.insert(make_pair(string("GC_26"),(((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))/pow(sw,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_26"]);
      m_local_complex_constants.insert(make_pair(string("GC_26"),(((pow(cw,2.0)*pow(ee,2.0))*complex(0.0,1.0))/pow(sw,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_26"]);
      m_local_complex_constants.insert(make_pair(string("GC_27"),((((-2.0*pow(cw,2.0))*pow(ee,2.0))*complex(0.0,1.0))/pow(sw,2.0))));
      DEBUG_VAR((m_local_complex_constants)["GC_27"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_3"),(-(ee*complex(0.0,1.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_3"]);
      m_local_complex_constants.insert(make_pair(string("GC_2"),(((2.0*ee)*complex(0.0,1.0))/3.0)));
      DEBUG_VAR((m_local_complex_constants)["GC_2"]);
      m_local_complex_constants.insert(make_pair(string("GC_2"),(((2.0*ee)*complex(0.0,1.0))/3.0)));
      DEBUG_VAR((m_local_complex_constants)["GC_2"]);
      m_local_complex_constants.insert(make_pair(string("GC_2"),(((2.0*ee)*complex(0.0,1.0))/3.0)));
      DEBUG_VAR((m_local_complex_constants)["GC_2"]);
      m_local_complex_constants.insert(make_pair(string("GC_1"),((-(ee*complex(0.0,1.0)))/3.0)));
      DEBUG_VAR((m_local_complex_constants)["GC_1"]);
      m_local_complex_constants.insert(make_pair(string("GC_1"),((-(ee*complex(0.0,1.0)))/3.0)));
      DEBUG_VAR((m_local_complex_constants)["GC_1"]);
      m_local_complex_constants.insert(make_pair(string("GC_1"),((-(ee*complex(0.0,1.0)))/3.0)));
      DEBUG_VAR((m_local_complex_constants)["GC_1"]);
      m_local_complex_constants.insert(make_pair(string("GC_12"),(complex(0.0,1.0)*G)));
      DEBUG_VAR((m_local_complex_constants)["GC_12"]);
      m_local_complex_constants.insert(make_pair(string("GC_12"),(complex(0.0,1.0)*G)));
      DEBUG_VAR((m_local_complex_constants)["GC_12"]);
      m_local_complex_constants.insert(make_pair(string("GC_12"),(complex(0.0,1.0)*G)));
      DEBUG_VAR((m_local_complex_constants)["GC_12"]);
      m_local_complex_constants.insert(make_pair(string("GC_12"),(complex(0.0,1.0)*G)));
      DEBUG_VAR((m_local_complex_constants)["GC_12"]);
      m_local_complex_constants.insert(make_pair(string("GC_12"),(complex(0.0,1.0)*G)));
      DEBUG_VAR((m_local_complex_constants)["GC_12"]);
      m_local_complex_constants.insert(make_pair(string("GC_12"),(complex(0.0,1.0)*G)));
      DEBUG_VAR((m_local_complex_constants)["GC_12"]);
      m_local_complex_constants.insert(make_pair(string("GC_33"),(((CKM1x1*ee)*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_33"]);
      m_local_complex_constants.insert(make_pair(string("GC_34"),(((CKM1x2*ee)*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_34"]);
      m_local_complex_constants.insert(make_pair(string("GC_35"),(((CKM2x1*ee)*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_35"]);
      m_local_complex_constants.insert(make_pair(string("GC_36"),(((CKM2x2*ee)*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_36"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_79"),(((ee*complex(0.0,1.0))*conj(CKM1x1))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_79"]);
      m_local_complex_constants.insert(make_pair(string("GC_81"),(((ee*complex(0.0,1.0))*conj(CKM2x1))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_81"]);
      m_local_complex_constants.insert(make_pair(string("GC_80"),(((ee*complex(0.0,1.0))*conj(CKM1x2))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_80"]);
      m_local_complex_constants.insert(make_pair(string("GC_82"),(((ee*complex(0.0,1.0))*conj(CKM2x2))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_82"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_32"),((ee*complex(0.0,1.0))/(sw*sqrt(2.0)))));
      DEBUG_VAR((m_local_complex_constants)["GC_32"]);
      m_local_complex_constants.insert(make_pair(string("GC_45"),((((-2.0*ee)*complex(0.0,1.0))*sw)/(3.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_45"]);
      m_local_complex_constants.insert(make_pair(string("GC_49"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(6.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_49"]);
      m_local_complex_constants.insert(make_pair(string("GC_45"),((((-2.0*ee)*complex(0.0,1.0))*sw)/(3.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_45"]);
      m_local_complex_constants.insert(make_pair(string("GC_49"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(6.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_49"]);
      m_local_complex_constants.insert(make_pair(string("GC_45"),((((-2.0*ee)*complex(0.0,1.0))*sw)/(3.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_45"]);
      m_local_complex_constants.insert(make_pair(string("GC_49"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(6.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_49"]);
      m_local_complex_constants.insert(make_pair(string("GC_44"),(((ee*complex(0.0,1.0))*sw)/(3.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_44"]);
      m_local_complex_constants.insert(make_pair(string("GC_48"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(6.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_48"]);
      m_local_complex_constants.insert(make_pair(string("GC_44"),(((ee*complex(0.0,1.0))*sw)/(3.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_44"]);
      m_local_complex_constants.insert(make_pair(string("GC_48"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(6.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_48"]);
      m_local_complex_constants.insert(make_pair(string("GC_44"),(((ee*complex(0.0,1.0))*sw)/(3.0*cw))));
      DEBUG_VAR((m_local_complex_constants)["GC_44"]);
      m_local_complex_constants.insert(make_pair(string("GC_48"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))-(((ee*complex(0.0,1.0))*sw)/(6.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_48"]);
      m_local_complex_constants.insert(make_pair(string("GC_52"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_52"]);
      m_local_complex_constants.insert(make_pair(string("GC_52"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_52"]);
      m_local_complex_constants.insert(make_pair(string("GC_52"),((((cw*ee)*complex(0.0,1.0))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_52"]);
      m_local_complex_constants.insert(make_pair(string("GC_46"),(((ee*complex(0.0,1.0))*sw)/cw)));
      DEBUG_VAR((m_local_complex_constants)["GC_46"]);
      m_local_complex_constants.insert(make_pair(string("GC_51"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_51"]);
      m_local_complex_constants.insert(make_pair(string("GC_46"),(((ee*complex(0.0,1.0))*sw)/cw)));
      DEBUG_VAR((m_local_complex_constants)["GC_46"]);
      m_local_complex_constants.insert(make_pair(string("GC_51"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_51"]);
      m_local_complex_constants.insert(make_pair(string("GC_46"),(((ee*complex(0.0,1.0))*sw)/cw)));
      DEBUG_VAR((m_local_complex_constants)["GC_46"]);
      m_local_complex_constants.insert(make_pair(string("GC_51"),(((-((cw*ee)*complex(0.0,1.0)))/(2.0*sw))+(((ee*complex(0.0,1.0))*sw)/(2.0*cw)))));
      DEBUG_VAR((m_local_complex_constants)["GC_51"]);
      msg_Debugging() << setprecision(6);
    }
    
    void vertices_0 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)15,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)15,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_46",m_local_complex_constants["GC_46"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_51",m_local_complex_constants["GC_51"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)13,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)13,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_46",m_local_complex_constants["GC_46"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_51",m_local_complex_constants["GC_51"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)11,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)11,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_46",m_local_complex_constants["GC_46"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_51",m_local_complex_constants["GC_51"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)16,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)16,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_52",m_local_complex_constants["GC_52"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)14,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)14,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_52",m_local_complex_constants["GC_52"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)12,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)12,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_52",m_local_complex_constants["GC_52"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_44",m_local_complex_constants["GC_44"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_48",m_local_complex_constants["GC_48"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_44",m_local_complex_constants["GC_44"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_48",m_local_complex_constants["GC_48"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_44",m_local_complex_constants["GC_44"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_48",m_local_complex_constants["GC_48"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_45",m_local_complex_constants["GC_45"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_49",m_local_complex_constants["GC_49"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
    }
    void vertices_1 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_45",m_local_complex_constants["GC_45"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_49",m_local_complex_constants["GC_49"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_45",m_local_complex_constants["GC_45"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_49",m_local_complex_constants["GC_49"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)16,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)15,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)14,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)13,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)12,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)11,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)15,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)16,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)13,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)14,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)11,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)12,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_82",m_local_complex_constants["GC_82"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
    }
    void vertices_2 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_80",m_local_complex_constants["GC_80"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_81",m_local_complex_constants["GC_81"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_79",m_local_complex_constants["GC_79"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_32",m_local_complex_constants["GC_32"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_36",m_local_complex_constants["GC_36"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_35",m_local_complex_constants["GC_35"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_34",m_local_complex_constants["GC_34"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_33",m_local_complex_constants["GC_33"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_12",m_local_complex_constants["GC_12"]) );
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_12",m_local_complex_constants["GC_12"]) );
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
    }
    void vertices_3 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_12",m_local_complex_constants["GC_12"]) );
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_12",m_local_complex_constants["GC_12"]) );
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_12",m_local_complex_constants["GC_12"]) );
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_12",m_local_complex_constants["GC_12"]) );
      m_v.back().Color.push_back(Color_Function(cf::T,3,2,1));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_1",m_local_complex_constants["GC_1"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)3,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_1",m_local_complex_constants["GC_1"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)1,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_1",m_local_complex_constants["GC_1"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_2",m_local_complex_constants["GC_2"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)4,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_2",m_local_complex_constants["GC_2"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)2,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_2",m_local_complex_constants["GC_2"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
    }
    void vertices_4 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)15,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)15,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_3",m_local_complex_constants["GC_3"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)13,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)13,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_3",m_local_complex_constants["GC_3"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)11,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)11,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_3",m_local_complex_constants["GC_3"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("FFV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_27",m_local_complex_constants["GC_27"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_26",m_local_complex_constants["GC_26"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_26",m_local_complex_constants["GC_26"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVVV4");
      m_v.back().Lorentz.push_back("VVVV3");
      m_v.back().Lorentz.push_back("VVVV2");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_72",m_local_complex_constants["GC_72"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_56",m_local_complex_constants["GC_56"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_55",m_local_complex_constants["GC_55"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_56",m_local_complex_constants["GC_56"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_42",m_local_complex_constants["GC_42"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_42",m_local_complex_constants["GC_42"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_43",m_local_complex_constants["GC_43"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVVV4");
      m_v.back().Lorentz.push_back("VVVV3");
      m_v.back().Lorentz.push_back("VVVV2");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_57",m_local_complex_constants["GC_57"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
    }
    void vertices_5 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_8",m_local_complex_constants["GC_8"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_9",m_local_complex_constants["GC_9"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_58",m_local_complex_constants["GC_58"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_10",m_local_complex_constants["GC_10"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_9",m_local_complex_constants["GC_9"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_50",m_local_complex_constants["GC_50"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_51",m_local_complex_constants["GC_51"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_53",m_local_complex_constants["GC_53"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_47",m_local_complex_constants["GC_47"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_54",m_local_complex_constants["GC_54"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_25",m_local_complex_constants["GC_25"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_24",m_local_complex_constants["GC_24"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_24",m_local_complex_constants["GC_24"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVVV4");
      m_v.back().Lorentz.push_back("VVVV3");
      m_v.back().Lorentz.push_back("VVVV2");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)23,0) );
      m_v.back().cpl.push_back( - ATOOLS::Kabbala("GC_37",m_local_complex_constants["GC_37"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
    }
    void vertices_6 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_6",m_local_complex_constants["GC_6"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_5",m_local_complex_constants["GC_5"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_5",m_local_complex_constants["GC_5"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVVV4");
      m_v.back().Lorentz.push_back("VVVV3");
      m_v.back().Lorentz.push_back("VVVV2");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_63",m_local_complex_constants["GC_63"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_23",m_local_complex_constants["GC_23"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_23",m_local_complex_constants["GC_23"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_23",m_local_complex_constants["GC_23"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_31",m_local_complex_constants["GC_31"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_28",m_local_complex_constants["GC_28"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_30",m_local_complex_constants["GC_30"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_29",m_local_complex_constants["GC_29"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_66",m_local_complex_constants["GC_66"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_41",m_local_complex_constants["GC_41"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_40",m_local_complex_constants["GC_40"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
    }
    void vertices_7 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,0) );
      m_v.back().cpl.push_back( - ATOOLS::Kabbala("GC_3",m_local_complex_constants["GC_3"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_31",m_local_complex_constants["GC_31"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_28",m_local_complex_constants["GC_28"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_29",m_local_complex_constants["GC_29"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_30",m_local_complex_constants["GC_30"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_65",m_local_complex_constants["GC_65"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_39",m_local_complex_constants["GC_39"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)24,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_40",m_local_complex_constants["GC_40"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_77",m_local_complex_constants["GC_77"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_76",m_local_complex_constants["GC_76"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_78",m_local_complex_constants["GC_78"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFSL");
      m_v.back().Lorentz.push_back("FFSR");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_19",m_local_complex_constants["GC_19"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_18",m_local_complex_constants["GC_18"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFSL");
      m_v.back().Lorentz.push_back("FFSR");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_74",m_local_complex_constants["GC_74"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
    }
    void vertices_8 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_75",m_local_complex_constants["GC_75"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_73",m_local_complex_constants["GC_73"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFSL");
      m_v.back().Lorentz.push_back("FFSR");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)6,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)5,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_17",m_local_complex_constants["GC_17"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_16",m_local_complex_constants["GC_16"]) );
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Color.push_back(Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFSL");
      m_v.back().Lorentz.push_back("FFSR");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_14",m_local_complex_constants["GC_14"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_15",m_local_complex_constants["GC_15"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_14",m_local_complex_constants["GC_14"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_15",m_local_complex_constants["GC_15"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_14",m_local_complex_constants["GC_14"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_15",m_local_complex_constants["GC_15"]) );
      m_v.back().Color.push_back(Color_Function(cf::F,-1,1,4,new Color_Function(cf::F,2,3,-1)));
      m_v.back().Color.push_back(Color_Function(cf::F,-1,1,4,new Color_Function(cf::F,2,3,-1)));
      m_v.back().Color.push_back(Color_Function(cf::F,-1,1,3,new Color_Function(cf::F,2,4,-1)));
      m_v.back().Color.push_back(Color_Function(cf::F,-1,1,3,new Color_Function(cf::F,2,4,-1)));
      m_v.back().Color.push_back(Color_Function(cf::F,-1,1,3,new Color_Function(cf::F,2,4,-1)));
      m_v.back().Color.push_back(Color_Function(cf::F,-1,1,2,new Color_Function(cf::F,3,4,-1)));
      m_v.back().Lorentz.push_back("VVVV4");
      m_v.back().Lorentz.push_back("VVVV3");
      m_v.back().Lorentz.push_back("VVVV4");
      m_v.back().Lorentz.push_back("VVVV2");
      m_v.back().Lorentz.push_back("VVVV3");
      m_v.back().Lorentz.push_back("VVVV2");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 2;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)21,0) );
      m_v.back().cpl.push_back( - ATOOLS::Kabbala("GC_13",m_local_complex_constants["GC_13"]) );
      m_v.back().Color.push_back(Color_Function(cf::F,1,2,3));
      m_v.back().Lorentz.push_back("VVV");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 1;
      m_v.back().order[1]    = 0;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_4",m_local_complex_constants["GC_4"]) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_3",m_local_complex_constants["GC_3"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VSS2");
      m_v.back().Lorentz.push_back("VSS1");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)22,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_7",m_local_complex_constants["GC_7"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("VVSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_60",m_local_complex_constants["GC_60"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_59",m_local_complex_constants["GC_59"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_59",m_local_complex_constants["GC_59"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 1;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_22",m_local_complex_constants["GC_22"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 2;
    }
    void vertices_9 () {
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_20",m_local_complex_constants["GC_20"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)25,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_20",m_local_complex_constants["GC_20"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_21",m_local_complex_constants["GC_21"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,1) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)42,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_20",m_local_complex_constants["GC_20"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().AddParticle( ATOOLS::Flavour((kf_code)41,0) );
      m_v.back().cpl.push_back( ATOOLS::Kabbala("GC_22",m_local_complex_constants["GC_22"]) );
      m_v.back().Color.push_back(Color_Function(cf::None));
      m_v.back().Lorentz.push_back("SSSS");
      m_v.back().order.resize(3);
      m_v.back().order[0]    = 0;
      m_v.back().order[1]    = 0;
      m_v.back().order[2]    = 2;
    }
    void InitVertices()
    {
      ParamInit(); 
      vertices_0();
      vertices_1();
      vertices_2();
      vertices_3();
      vertices_4();
      vertices_5();
      vertices_6();
      vertices_7();
      vertices_8();
      vertices_9();
    }

    void ResetVerticesWithEWParameters(const EWParameters& params)
    {
      ClearInteractionModel();
      /// Set parameters to their run value before re-initing the vertices
      (*p_complexconstants)[std::string("csin2_thetaW")] = params.m_sw2_r;
      (*p_complexconstants)[std::string("ccos2_thetaW")] = params.m_cw2_r;
      (*p_complexconstants)[std::string("cvev")]         = params.m_cvev_r;
      (*p_constants)[std::string("mZ")]                  = params.m_mz_r;
      (*p_constants)[std::string("mW+")]                 = params.m_mw_r;
      (*p_constants)[std::string("mh0")]                 = params.m_mh0_r;
      (*p_constants)[std::string("mt")]                  = params.m_mt_r;
      (*p_constants)[std::string("alpha_QED")]           = params.m_aew_r;
      InitializeInteractionModel();
    }

    void ClearInteractionModel()
    {
      m_v.clear();
      m_ov.clear();
      m_fls.clear();
      m_vmap.clear();
      m_vtable.clear();
    }

    size_t IndexOfOrderKey(const std::string& key) const
    {
      PRINT_VAR(key);
      if(key == "SMGold")
	return 2;
      else return Model_Base::IndexOfOrderKey(key);
    }

    
  };
  
}

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
