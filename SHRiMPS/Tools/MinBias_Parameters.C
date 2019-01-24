#include "SHRiMPS/Tools/MinBias_Parameters.H" 
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHRIMPS;

MinBias_Parameters SHRIMPS::MBpars;

MinBias_Parameters::MinBias_Parameters()
{
  p_ffs      = new std::list<Form_Factor *>;
  p_eikonals = new std::list<Omega_ik *>;
  p_xsecs    = new XSecs_Container();
}

MinBias_Parameters::~MinBias_Parameters() {
  Reset();
}

void MinBias_Parameters::Reset() {
  while (!p_eikonals->empty()) {
    delete p_eikonals->back();
    p_eikonals->pop_back();
  }
  p_eikonals->clear();
  delete p_eikonals;
  while (!p_ffs->empty()) {
    delete p_ffs->back();
    p_ffs->pop_back();
  }
  p_ffs->clear();
  delete p_ffs;
  delete p_xsecs;
}

void MinBias_Parameters::SetXSecs(Cross_Sections * xsecs) {
  p_xsecs->xs_tot = xsecs->SigmaTot();
  p_xsecs->xs_in  = xsecs->SigmaInel();
  p_xsecs->xs_el  = xsecs->SigmaEl();
  p_xsecs->xs_SD  = xsecs->SigmaSD();
  p_xsecs->xs_DD  = xsecs->SigmaDD();
}


void MinBias_Parameters::Init() {
  RegisterDefaults();
  FillRunParameters();
  FillFormFactorParameters();
  FillEikonalParameters();
  FillLadderParameters();
  FillShowerLinkParameters();
}

void MinBias_Parameters::RegisterDefaults() const
{
  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  s["Shrimps_Mode"].SetDefault("Inelastic");
  s["MB_Weight_Mode"].SetDefault("Unweighted");
  s["bmax"].SetDefault(10.0);
  s["accu"].SetDefault(5e-4);
  s["FF_Form"].SetDefault("dipole");
  s["Lambda2"].SetDefault(1.7);
  s["beta02(mb)"].SetDefault(25.0);
  s["kappa"].SetDefault(0.6);
  s["xi"].SetDefault(0.2);
  s["bsteps_FF"].SetDefault(64);
  s["Absorption"].SetDefault("factorial");
  s["deltaY"].SetDefault(1.5);
  s["lambda"].SetDefault(0.5);
  s["Delta"].SetDefault(0.3);
}

void MinBias_Parameters::FillRunParameters() {
  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  std::string runmode = s["Shrimps_Mode"].Get<std::string>();
  if (runmode==std::string("TestShrimps") || runmode==std::string("Test")) 
    m_runmode = m_run_params.runmode = run_mode::test;
  if (runmode==std::string("Xsecs") || runmode==std::string("XSecs")) 
    m_runmode = m_run_params.runmode = run_mode::xsecs_only;
  else if (runmode==std::string("Elastic")) 
    m_runmode = m_run_params.runmode = run_mode::elastic_events;
  else if (runmode==std::string("Single-diffractive")) 
    m_runmode = m_run_params.runmode = run_mode::single_diffractive_events;
  else if (runmode==std::string("Double-diffractive")) 
    m_runmode = m_run_params.runmode = run_mode::double_diffractive_events;
  else if (runmode==std::string("Quasi-elastic")) 
    m_runmode = m_run_params.runmode = run_mode::quasi_elastic_events;
  else if (runmode==std::string("Inelastic")) 
    m_runmode = m_run_params.runmode = run_mode::inelastic_events;
  else if (runmode==std::string("All")) 
    m_runmode = m_run_params.runmode = run_mode::all_min_bias;
  else if (runmode==std::string("Underlying")) 
    m_runmode = m_run_params.runmode = run_mode::underlying_event;

  msg_Out()<<METHOD<<"(mode = "<<runmode<<" -> "<<int(m_runmode)<<").\n";
  
  std::string weightmode = s["MB_Weight_Mode"].Get<std::string>();
  if (weightmode==std::string("Unweighted")) 
    m_run_params.weightmode = weight_mode::unweighted;
  else if (weightmode==std::string("Weighted")) 
    m_run_params.weightmode = weight_mode::weighted;

  m_originalY = log(ATOOLS::rpa->gen.Ecms()/
		    ATOOLS::Flavour(kf_p_plus).HadMass());
  m_NGWstates = (m_runmode!=run_mode::test)?2:1;
  m_bmax      = s["bmax"].Get<double>();
  m_accu      = s["accu"].Get<double>();
}

void MinBias_Parameters::FillFormFactorParameters() {
  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  std::string ffform = s["FF_Form"].Get<std::string>();
  if (ffform==std::string("dipole") && m_runmode!=run_mode::test) 
    m_ff_params.form = ff_form::dipole;
  else                               
    m_ff_params.form = ff_form::Gauss;
  m_ff_params.norm        = 1./sqrt(m_NGWstates);
  m_ff_params.Lambda2     = s["Lambda2"].Get<double>();
  m_ff_params.beta02      = 
    sqrt(1.e9*s["beta02(mb)"].Get<double>()/ATOOLS::rpa->Picobarn());
  // kappa will obtain a sign for the second GW state - this is hardwired
  // in the initialization routine in Shrimps.  
  m_ff_params.kappa       = s["kappa"].Get<double>();
  m_ff_params.xi          = s["xi"].Get<double>();
  m_ff_params.bmax        = m_bmax;
  m_ff_params.accu        = m_accu;
  m_ff_params.bsteps      = s["bsteps_FF"].Get<int>();
}

void MinBias_Parameters::FillEikonalParameters() {
  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  std::string absorption = s["Absorption"].Get<std::string>();
  if (absorption==std::string("exponential"))
    m_eik_params.absorp  = absorption::exponential;
  else
    m_eik_params.absorp  = absorption::factorial;
  m_eik_params.originalY = m_originalY;
  m_eik_params.cutoffY   = s["deltaY"].Get<double>();
  m_eik_params.Ymax      = m_eik_params.originalY - m_eik_params.cutoffY; 
  m_eik_params.lambda    = (m_runmode!=run_mode::test)?
    s["lambda"].Get<double>():0.;
  m_eik_params.Delta     = s["Delta"].Get<double>();
  m_eik_params.beta02    = m_ff_params.beta02;
  m_eik_params.bmax      = 2.*m_bmax;
  m_eik_params.accu      = m_accu;
}

void MinBias_Parameters::UpdateForNewEnergy(const double & energy) {
  m_originalY = log(energy/ATOOLS::Flavour(kf_p_plus).HadMass());
  m_eik_params.originalY = m_originalY;
  m_eik_params.Ymax      = m_eik_params.originalY - m_eik_params.cutoffY; 
}


void MinBias_Parameters::FillLadderParameters() {
}

void MinBias_Parameters::FillShowerLinkParameters() {
}


