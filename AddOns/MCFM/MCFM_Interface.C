#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MCFM/CXX_Interface.h"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include <ctype.h>

using namespace PHASIC; 
using namespace ATOOLS;

namespace SHERPA {

  inline std::string str_tolower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
		   [](unsigned char c){ return tolower(c); });
    return s;
  }

  class MCFM_Interface: public PHASIC::ME_Generator_Base {
  private:

    static MCFM::CXX_Interface s_mcfm;
    static MODEL::Running_AlphaS *p_as;

  public:

    MCFM_Interface(): ME_Generator_Base("MCFM")
    {
      rpa->gen.AddCitation(1,s_mcfm.GetReferences());
      msg_Info()<<s_mcfm.GetStartupMessage();
    }

    ~MCFM_Interface()
    {
      msg_Info()<<s_mcfm.GetFinishMessage(1);
    }

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr)
    {
      DEBUG_FUNC("");
      p_as=(MODEL::Running_AlphaS*)model->GetScalarFunction("alpha_S");
      std::string pdname(rpa->gen.Variable("SHERPA_CPP_PATH")+"/process.DAT");
      if (!FileExists(pdname))
	Copy(MCFM_PATH+std::string("/share/MCFM/process.DAT"),pdname);
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
      params["charm_yukawa"]=ToString(std::abs(model->ComplexConstant("yuk(4)")),16);
      params["bottom_yukawa"]=ToString(std::abs(model->ComplexConstant("yuk(5)")),16);
      params["top_yukawa"]=ToString(std::abs(model->ComplexConstant("yuk(6)")),16);
      params["tau_yukawa"]=ToString(std::abs(model->ComplexConstant("yuk(15)")),16);
      params["ew_scheme"]="5";
      params["alpha_EM"]=ToString(model->ScalarConstant("alpha_QED"),16);
      params["Gf"]=ToString(1.0/sqrt(2.0)/std::abs(sqr(model->ComplexConstant("cvev"))),16);
      params["sin2_thetaW"]=ToString(std::abs(model->ComplexConstant("csin2_thetaW")),16);
      params["CKM_u_d"]=ToString(model->ComplexConstant("CKM_0_0"),16);
      params["CKM_u_s"]=ToString(model->ComplexConstant("CKM_1_0"),16);
      params["CKM_u_b"]=ToString(model->ComplexConstant("CKM_2_0"),16);
      params["CKM_c_d"]=ToString(model->ComplexConstant("CKM_0_1"),16);
      params["CKM_c_s"]=ToString(model->ComplexConstant("CKM_1_1"),16);
      params["CKM_c_b"]=ToString(model->ComplexConstant("CKM_2_1"),16);
      params["order_alpha_S"]=ToString(MODEL::as->Order()+1);
      params["alpha_S"]=ToString(model->ScalarConstant("alpha_S"),16);
      params["scale"]=ToString(Flavour(kf_Z).Mass(),16);
      s_mcfm.SetVerbose(msg->LevelIsDebugging());
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

    inline static void SetAlpha(const double &as,const double &aqed)
    {
      s_mcfm.SetAlphaS(as);
      // s_mcfm.SetAlphaQED(aqed);
    }

  }; // end of class MCFM_Interface

  class MCFM_Virtual: public PHASIC::Virtual_ME2_Base {
  private:

    MCFM::Process *p_proc;
    std::vector<MCFM::FourVec> m_p;

  public:

    MCFM_Virtual(const PHASIC::Process_Info& pi,
		 const ATOOLS::Flavour_Vector& flavs,int pid):
      Virtual_ME2_Base(pi,flavs),
      p_proc(MCFM_Interface::GetMCFM().GetProcess(pid))
    {
      rpa->gen.AddCitation(1,p_proc->GetReferences());
      m_p.resize(flavs.size());
      m_mode=1;
      m_drmode=p_proc->GetScheme();
    }

    void SetPoleCheck(const int check)
    {
      p_proc->SetPoleCheck(check);
    }

    void Calc(const ATOOLS::Vec4D_Vector &p)
    {
      MCFM_Interface::GetMCFM().
	SetVerbose(msg->LevelIsDebugging());
      for (size_t i(0);i<p.size();++i)
	for (size_t j(0);j<4;++j) m_p[i][j]=p[i][j];
      double ason2pi(MCFM_Interface::SetMuR2(m_mur2)/(2.*M_PI));
      p_proc->Calc(m_p,1);
      const std::vector<double> &res(p_proc->GetResult());
      m_res.Finite()=res[0]/ason2pi;
      m_res.IR()=res[1]/ason2pi;
      m_res.IR2()=res[2]/ason2pi;
      m_born=res[3]/p_proc->GetSymmetryFactor();
    }

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.*M_PI;// MSbar scheme
    }

  };// end of class MCFM_Virtual

  class MCFM_Born: public PHASIC::Tree_ME2_Base {

    int m_order_ew, m_order_qcd;
    MCFM::Process *p_proc;
    std::vector<MCFM::FourVec> m_p;

  public:

    MCFM_Born(const PHASIC::External_ME_Args &args,const int &pid):
      Tree_ME2_Base(args),
      p_proc(MCFM_Interface::GetMCFM().GetProcess(pid)) {
      rpa->gen.AddCitation(1,p_proc->GetReferences());
      m_p.resize(args.Flavours().size());
      m_order_qcd=args.m_orders[0];
      m_order_ew=args.m_orders[1];
    }

    double Calc(const ATOOLS::Vec4D_Vector &p)
    {
      MCFM_Interface::GetMCFM().
	SetVerbose(msg->LevelIsDebugging());
      for (size_t i(0);i<p.size();++i)
	for (size_t j(0);j<4;++j) m_p[i][j]=p[i][j];
      MCFM_Interface::SetAlpha(AlphaQCD(),AlphaQED());
      p_proc->Calc(m_p,0);
      const std::vector<double> &res(p_proc->GetResult());
      return res[3];
    }

    int OrderQCD(const int &id=-1) const { return m_order_qcd; }
    int OrderEW(const int &id=-1) const  { return m_order_ew;  }

  };// end of class MCFM_Born

} // end of namespace MCFM

using namespace SHERPA;

MCFM::CXX_Interface MCFM_Interface::s_mcfm(0);
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
  if (!(pi.m_fi.m_nlotype&nlo_type::loop)) return NULL;
  if (pi.m_fi.m_nlocpl[1]!=0.) return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  std::vector<int> ids(fl.size());
  for (size_t i(0);i<fl.size();++i) ids[i]=(long int)(fl[i]);
  MCFM::Process_Info mpi(ids,pi.m_ii.m_ps.size(),
			 pi.m_maxcpl[0],pi.m_maxcpl[1]+addcpl);
  std::string modelname(str_tolower(MODEL::s_model->Name()));
  if (modelname=="smehc") modelname="heft";
  mpi.m_model=str_tolower(modelname);
  DecayInfo_Vector decins(pi.m_fi.GetDecayInfos());
  for (size_t i(0);i<decins.size();++i) {
    mpi.m_decids.push_back(decins[i]->m_id);
    mpi.m_decfls.push_back((long int)(decins[i]->m_fl));
  }
  int pid(MCFM_Interface::GetMCFM().InitializeProcess(mpi));
  if (pid>=0) return new MCFM_Virtual(pi,fl,pid);
  return NULL;
}

DECLARE_TREEME2_GETTER(MCFM_Born,"MCFM_Born")

Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,External_ME_Args,MCFM_Born>::
operator()(const External_ME_Args &args) const
{
  if (args.m_source!="MCFM") return NULL;
  Flavour_Vector fl(args.Flavours());
  std::vector<int> ids(fl.size());
  for (size_t i(0);i<fl.size();++i) ids[i]=(long int)(fl[i]);
  MCFM::Process_Info mpi(ids,args.m_inflavs.size(),
			 args.m_orders[0],args.m_orders[1]);
  int pid(MCFM_Interface::GetMCFM().InitializeProcess(mpi));
  if (pid>=0) return new MCFM_Born(args,pid);
  return NULL;
}
