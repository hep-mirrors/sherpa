#include "SHRiMPS/Tools/MinBias_Parameters.H" 
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace SHRIMPS;
using namespace ATOOLS;

MinBias_Parameters SHRIMPS::MBpars;

MinBias_Parameters::MinBias_Parameters() { }

MinBias_Parameters::~MinBias_Parameters() {
  Reset();
}

void MinBias_Parameters::Reset() {
  ResetEikonals(0);
  while (!m_ffs.empty()) {
    delete m_ffs.back();
    m_ffs.pop_back();
  }
  m_ffs.clear();
}

void MinBias_Parameters::ResetEikonals(const size_t size) {
  if (!m_eikonals.empty()) {
    while (!m_eikonals.empty()) {
      while (!m_eikonals.back().empty()) {
	delete m_eikonals.back().back();
	m_eikonals.back().pop_back();
      }
      m_eikonals.pop_back();      
    }
    m_eikonals.clear();
  }
  m_eikonals.resize(size);
  for (size_t i=0;i<size;i++) m_eikonals[i].resize(size);
}


void MinBias_Parameters::SetXSecs(Cross_Sections * xsecs) {
  m_xsecs.xs_tot   = xsecs->SigmaTot();
  m_xsecs.xs_in    = xsecs->SigmaInel();
  m_xsecs.xs_el    = xsecs->SigmaEl();
  m_xsecs.xs_SD[0] = xsecs->SigmaSD(0);
  m_xsecs.xs_SD[1] = xsecs->SigmaSD(1);
  m_xsecs.xs_DD    = xsecs->SigmaDD();
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
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  s["Mode"].SetDefault("Inelastic");
  s["MB_Weight_Mode"].SetDefault("Unweighted");
  s["bmax"].SetDefault(10.0);
  s["accu"].SetDefault(5e-4);
  s["bsteps_FF"].SetDefault(64);
  s["Absorption"].SetDefault("exponential");
  s["FF_Form"].SetDefault("dipole");
  s["deltaY"].SetDefault(0.1);
  //s["deltaY"].SetDefault(0.05);
  s["beta02(mb)"].SetDefault(20.0);
  s["Lambda2"].SetDefault(1.2);
  s["kappa"].SetDefault(0.55);
  s["xi"].SetDefault(0.2);
  s["lambda"].SetDefault(0.25);
  s["Delta"].SetDefault(0.4);
  s["Q_0^2"].SetDefault(1.0);
  s["Q_as^2"].SetDefault(1.0);
  s["KT_shower"].SetDefault(5.);
  s["Collinear_KT_min"].SetDefault(1.);
  s["Incl_Tune"].SetDefault("None");
}

void MinBias_Parameters::FillRunParameters() {
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  std::string runmode = s["Mode"].Get<std::string>();
  if (runmode==std::string("TestShrimps") || runmode==std::string("Test")) 
    m_runmode = m_run_params.runmode = run_mode::test;
  else if (runmode==std::string("Xsecs") || runmode==std::string("XSecs")) 
    m_runmode = m_run_params.runmode = run_mode::xsecs_only;
  else if (runmode==std::string("Elastic")) 
    m_runmode = m_run_params.runmode = run_mode::elastic_events;
  else if (runmode==std::string("Soft-Diffractive")) 
    m_runmode = m_run_params.runmode = run_mode::soft_diffractive_events;
  else if (runmode==std::string("Quasi-elastic")) 
    m_runmode = m_run_params.runmode = run_mode::quasi_elastic_events;
  else if (runmode==std::string("Inelastic")) 
    m_runmode = m_run_params.runmode = run_mode::inelastic_events;
  else if (runmode==std::string("DD")) 
    m_runmode = m_run_params.runmode = run_mode::double_diffractive_events;
  else if (runmode==std::string("SD")) 
    m_runmode = m_run_params.runmode = run_mode::single_diffractive_events;
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

  m_originalY = log(rpa->gen.Ecms()/
		    Flavour(kf_p_plus).HadMass());
  m_NGWstates = (m_runmode!=run_mode::test)?2:1;
  m_bmax      = s["bmax"].Get<double>();
  m_accu      = s["accu"].Get<double>();
}

void MinBias_Parameters::FillFormFactorParameters() {
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  std::string ffform = s["FF_Form"].Get<std::string>();
  if (ffform==std::string("dipole") && m_runmode!=run_mode::test) 
    m_ff_params.form = ff_form::dipole;
  else                               
    m_ff_params.form = ff_form::Gauss;
  m_ff_params.norm        = 1./sqrt(m_NGWstates);
  m_ff_params.Lambda2     = s["Lambda2"].Get<double>();
  m_ff_params.beta02      = 
    sqrt(1.e9*s["beta02(mb)"].Get<double>()/rpa->Picobarn());
  // kappa will obtain a sign for the second GW state - this is hardwired
  // in the initialization routine in Shrimps.  
  m_ff_params.kappa       = s["kappa"].Get<double>();
  m_ff_params.xi          = s["xi"].Get<double>();
  m_ff_params.bmax        = m_bmax;
  m_ff_params.accu        = m_accu;
  m_ff_params.bsteps      = s["bsteps_FF"].Get<int>();
  std::string incltune = s["Incl_Tune"].Get<std::string>();
  if (incltune==std::string("tune1")) {
      m_ff_params.beta02  = sqrt(1.e9*29.13/rpa->Picobarn());
      m_ff_params.Lambda2 = 1.683;
      m_ff_params.kappa   = 0.5966;
      m_ff_params.xi      = 0.1014;
  }
  if (incltune==std::string("tune2")) {
      m_ff_params.beta02  = sqrt(1.e9*15.02/rpa->Picobarn());
      m_ff_params.Lambda2 = 1.205;
      m_ff_params.kappa   = 0.4583;
      m_ff_params.xi      = 0.1992;
  }
  if (incltune!=std::string("tune1") && incltune!=std::string("tune1")) {
      msg_Out()<<METHOD<<": Unrecognised inclusive tune: "<<incltune<<", will fall back to default settings.\n";
  }
  else {
      msg_Out()<<METHOD<<": Using inclusive tune: "<<incltune<<"\n";
  }
  msg_Out()<<"    beta02(mb) = "<<sqr(m_ff_params.beta02)*rpa->Picobarn()/1.e9<<"\n";
  msg_Out()<<"    Lambda2    = "<<m_ff_params.Lambda2<<"\n";
  msg_Out()<<"    kappa      = "<<m_ff_params.kappa<<"\n";
  msg_Out()<<"    xi         = "<<m_ff_params.xi<<std::endl;
}

void MinBias_Parameters::FillEikonalParameters() {
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  std::string absorption = s["Absorption"].Get<std::string>();
  if (absorption==std::string("exponential"))
    m_eik_params.absorp  = absorption::exponential;
  else
    m_eik_params.absorp  = absorption::factorial;
  m_eik_params.originalY = m_originalY;
  //m_eik_params.cutoffY   = s["deltaY"].Get<double>();
  m_eik_params.cutoffY   = s["deltaY"].Get<double>()*m_originalY;
  m_eik_params.Ymax      = m_eik_params.originalY - m_eik_params.cutoffY;
  m_eik_params.lambda    = (m_runmode!=run_mode::test)?
    s["lambda"].Get<double>():0.;
  m_eik_params.Delta     = s["Delta"].Get<double>();
  m_eik_params.beta02    = m_ff_params.beta02;
  m_eik_params.bmax      = 2.*m_bmax;
  m_eik_params.accu      = m_accu;
  std::string incltune = s["Incl_Tune"].Get<std::string>();
  if (incltune==std::string("tune1")) {
      m_eik_params.absorp  = absorption::exponential;
      m_eik_params.cutoffY = 0.001358*m_originalY;
      m_eik_params.lambda  = 0.1782;
      m_eik_params.Delta   = 0.4988;
  }
  if (incltune==std::string("tune2")) {
      m_eik_params.absorp  = absorption::exponential;
      m_eik_params.cutoffY = 0.03272*m_originalY;
      m_eik_params.lambda  = 0.2566;
      m_eik_params.Delta   = 0.4403;
  }
  if (incltune!=std::string("tune1") && incltune!=std::string("tune1")) {
      msg_Out()<<METHOD<<": Unrecognised inclusive tune: "<<incltune<<", will fall back to default settings.\n";
  }
  else {
      msg_Out()<<METHOD<<": Using inclusive tune: "<<incltune<<"\n";
  }
  msg_Out()<<"    Absorption = "<<(m_eik_params.absorp==absorption::exponential?"exponential":"factorial")<<"\n";
  msg_Out()<<"    deltaY     = "<<m_eik_params.cutoffY/m_originalY<<"\n";
  msg_Out()<<"    lambda     = "<<m_eik_params.lambda<<"\n";
  msg_Out()<<"    Delta      = "<<m_eik_params.Delta<<std::endl;
}

void MinBias_Parameters::UpdateForNewEnergy(const double & energy) {
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  m_originalY = log(energy/Flavour(kf_p_plus).HadMass());
  m_eik_params.originalY = m_originalY;
  m_eik_params.cutoffY   = s["deltaY"].Get<double>()*m_originalY;
  m_eik_params.Ymax      = m_eik_params.originalY - m_eik_params.cutoffY;
}


void MinBias_Parameters::FillLadderParameters() {
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  m_ladder_params.Q02  = s["Q_0^2"].Get<double>();
  m_ladder_params.Qas2 = s["Q_as^2"].Get<double>();
}

void MinBias_Parameters::FillShowerLinkParameters() {
  const Scoped_Settings & s = Settings::GetMainSettings()["SHRIMPS"];
  m_showerlink_params.KT2min = sqr(s["KT_shower"].Get<double>());
  m_showerlink_params.CEKT2min = sqr(s["Collinear_KT_min"].Get<double>());
}


