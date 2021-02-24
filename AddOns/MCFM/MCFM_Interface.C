#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MCFM/CXX_Interface.h"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC; 
using namespace ATOOLS;

namespace SHERPA {

  class MCFM_Interface: public PHASIC::ME_Generator_Base {
  private:

    static MCFM::CXX_Interface s_mcfm;
    static MODEL::Running_AlphaS *p_as;

  public:

    MCFM_Interface(): ME_Generator_Base("MCFM") {}

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr)
    {
      DEBUG_FUNC("");
      p_as=(MODEL::Running_AlphaS*)model->GetScalarFunction("alpha_S");
      std::string pdname(rpa->gen.Variable("SHERPA_CPP_PATH")+"/process.DAT");
      if (!FileExists(pdname))
	Copy(MCFM_PATH+std::string("/share/process.DAT"),pdname);
      std::map<std::string,std::string> params;
      params["n_flav"]=ToString(Flavour(kf_jet).Size()/2,16);
      params["down_mass"]=ToString(Flavour(kf_d).Mass(),16);
      params["up_mass"]=ToString(Flavour(kf_u).Mass(),16);
      params["strange_mass"]=ToString(Flavour(kf_s).Mass(),16);
      params["charm_mass"]=ToString(Flavour(kf_c).Mass(),16);
      params["bottom_mass"]=ToString(Flavour(kf_b).Mass(),16);
      params["top_mass"]=ToString(Flavour(kf_t).Mass(),16);
      params["top_width"]=ToString(Flavour(kf_t).Width(),16);
      params["electron_mass"]=ToString(Flavour(kf_e).Mass(),16);
      params["muon_mass"]=ToString(Flavour(kf_mu).Mass(),16);
      params["tau_mass"]=ToString(Flavour(kf_tau).Mass(),16);
      params["tau_width"]=ToString(Flavour(kf_tau).Width(),16);
      params["H_mass"]=ToString(Flavour(kf_h0).Mass(),16);
      params["H_width"]=ToString(Flavour(kf_h0).Width(),16);
      params["Z_mass"]=ToString(Flavour(kf_Z).Mass(),16);
      params["Z_width"]=ToString(Flavour(kf_Z).Width(),16);
      params["W_mass"]=ToString(Flavour(kf_Wplus).Mass(),16);
      params["W_width"]=ToString(Flavour(kf_Wplus).Width(),16);
      params["charm_mass_square"]=ToString(sqr(Flavour(kf_c).Mass(true)),16);
      params["bottom_mass_square"]=ToString(sqr(Flavour(kf_b).Mass(true)),16);
      params["tau_mass_square"]=ToString(sqr(Flavour(kf_tau).Mass(true)),16);
      params["alpha_EM"]=ToString(model->ScalarConstant("alpha_QED"),16);
      params["Gf"]=ToString(1.0/sqrt(2.0)/std::abs(sqr(model->ComplexConstant("cvev"))),16);
      params["sin2_thetaW"]=ToString(std::abs(model->ComplexConstant("csin2_thetaW")),16);
      params["CKM_u_d"]=ToString(model->ComplexConstant("CKM_0_0").real(),16);
      params["CKM_u_s"]=ToString(model->ComplexConstant("CKM_0_1").real(),16);
      params["CKM_u_b"]=ToString(model->ComplexConstant("CKM_0_2").real(),16);
      params["CKM_c_d"]=ToString(model->ComplexConstant("CKM_1_0").real(),16);
      params["CKM_c_s"]=ToString(model->ComplexConstant("CKM_1_1").real(),16);
      params["CKM_c_b"]=ToString(model->ComplexConstant("CKM_1_2").real(),16);
      params["order_alpha_S"]=ToString(MODEL::as->Order()+1);
      params["alpha_S"]=ToString(model->ScalarConstant("alpha_S"),16);
      s_mcfm.Initialize(params);
      return true;
    }

    PHASIC::Process_Base *InitializeProcess
    (const PHASIC::Process_Info &pi, bool add) { return NULL; }

    int  PerformTests() { return 1; }
    bool NewLibraries() { return false; }

    inline static MCFM::CXX_Interface &GetMCFM() { return s_mcfm; }

    inline static double SetMuR2(const double &mur2)
    {
      double as((*p_as)(mur2));
      s_mcfm.SetMuR2(mur2);
      s_mcfm.SetAlphaS(as);
      return as;
    }

  }; // end of class MCFM_Interface

  class MCFM_Virtual: public PHASIC::Virtual_ME2_Base {
  private:

    int m_pid;
    std::vector<MCFM::FourVec> m_p;

  public:

    MCFM_Virtual(const PHASIC::Process_Info& pi,
		 const ATOOLS::Flavour_Vector& flavs,int pid):
      Virtual_ME2_Base(pi,flavs), m_pid(pid)
    {
      rpa->gen.AddCitation
	(1,"NLO matrix elements from MCFM \\cite{}.");
      m_p.resize(flavs.size());
      m_mode=1;
      m_drmode=MCFM_Interface::GetMCFM().GetScheme(pid);
    }

    void Calc(const ATOOLS::Vec4D_Vector &p)
    {
      for (size_t i(0);i<p.size();++i)
	for (size_t j(0);j<4;++j) m_p[i][j]=p[i][j];
      double ason2pi(MCFM_Interface::SetMuR2(m_mur2)/(2.*M_PI));
      MCFM_Interface::GetMCFM().Calc(m_pid,m_p,1);
      const std::vector<double> &res
	(MCFM_Interface::GetMCFM().GetResult(m_pid));
      m_res.Finite()=res[0]/ason2pi;
      m_res.IR()=res[1]/ason2pi;
      m_res.IR2()=res[2]/ason2pi;
      m_born=res[3];
    }

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.*M_PI;// MSbar scheme
    }

  };// end of class MCFM_Virtual

} // end of namespace MCFM

using namespace SHERPA;

MCFM::CXX_Interface MCFM_Interface::s_mcfm;
MODEL::Running_AlphaS *MCFM_Interface::p_as(NULL);

DECLARE_GETTER(MCFM_Interface,"MCFM",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,MCFM_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new MCFM_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,MCFM_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Interface to the MCFM loop ME generator"; 
}

DECLARE_VIRTUALME2_GETTER(MCFM_Virtual,"MCFM_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,MCFM_Virtual>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (MODEL::s_model->Name()!="SM") return NULL;
  if (!(pi.m_fi.m_nlotype&nlo_type::loop)) return NULL;
  if (pi.m_fi.m_nlocpl[1]!=0.) return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  std::vector<int> ids(fl.size());
  for (size_t i(0);i<fl.size();++i) ids[i]=(long int)(fl[i]);
  int pid(MCFM_Interface::GetMCFM().InitializeProcess(ids));
  if (pid>=0) return new MCFM_Virtual(pi,fl,pid);
  return NULL;
}
