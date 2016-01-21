#include "SHRiMPS/Tools/MinBias_Parameters.H" 
#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;

MinBias_Parameters SHRIMPS::MBpars;

MinBias_Parameters::MinBias_Parameters()
{
  p_ffs      = new std::list<Form_Factor *>;
  p_eikonals = new std::list<Omega_ik *>;
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
}

void MinBias_Parameters::Init(ATOOLS::Data_Reader * dr) {
  FillRunParameters(dr);
  FillFormFactorParameters(dr);
  FillEikonalParameters(dr);
  FillLadderParameters(dr);
  FillShowerLinkParameters(dr);
}

void MinBias_Parameters::FillRunParameters(ATOOLS::Data_Reader * dr) {
  std::string runmode =
    dr->GetValue<std::string>("Shrimps_Mode",std::string("Inelastic"));
  if (runmode==std::string("TestShrimps") || runmode==std::string("Test")) 
    m_runmode = m_run_params.runmode = run_mode::test;
  if (runmode==std::string("Xsecs")) 
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
  
  std::string weightmode =      
    dr->GetValue<std::string>("MB_Weight_Mode",std::string("Unweighted"));
  if (weightmode==std::string("Unweighted")) 
    m_run_params.weightmode = weight_mode::unweighted;
  else if (weightmode==std::string("Weighted")) 
    m_run_params.weightmode = weight_mode::weighted;

  m_originalY = log(ATOOLS::rpa->gen.Ecms()/
		    ATOOLS::Flavour(kf_p_plus).HadMass());
  m_NGWstates = (m_runmode!=run_mode::test)?2:1;
  m_bmax      = dr->GetValue<double>("bmax",10.);
  m_accu      = dr->GetValue<double>("accu",5.e-4);
}

void MinBias_Parameters::FillFormFactorParameters(ATOOLS::Data_Reader * dr) {
  std::string ffform = 
    dr->GetValue<std::string>("FF_Form",std::string("dipole"));
  if (ffform==std::string("dipole") && m_runmode!=run_mode::test) 
    m_ff_params.form = ff_form::dipole;
  else                               
    m_ff_params.form = ff_form::Gauss;
  m_ff_params.norm        = 1./sqrt(m_NGWstates);
  m_ff_params.Lambda2     = dr->GetValue<double>("Lambda2",1.7);
  m_ff_params.beta02      = 
    sqrt(1.e9*dr->GetValue<double>("beta02(mb)",20.0)/ATOOLS::rpa->Picobarn());
  // kappa will obtain a sign for the second GW state - this is hardwired
  // in the initialization routine in Shrimps.  
  m_ff_params.kappa       = dr->GetValue<double>("kappa",0.6);
  m_ff_params.xi          = dr->GetValue<double>("xi",0.2);
  m_ff_params.bmax        = m_bmax;
  m_ff_params.accu        = m_accu;
  m_ff_params.bsteps      = dr->GetValue<int>("bsteps_FF",64);
}

void MinBias_Parameters::FillEikonalParameters(ATOOLS::Data_Reader * dr) {
  std::string absorption = 
    dr->GetValue<std::string>("Absorption",std::string("factorial"));
  if (absorption==std::string("exponential"))
    m_eik_params.absorp  = absorption::exponential;
  else
    m_eik_params.absorp  = absorption::factorial;
  m_eik_params.originalY = m_originalY;
  m_eik_params.cutoffY   = dr->GetValue<double>("deltaY",1.5);
  m_eik_params.Ymax      = m_eik_params.originalY - m_eik_params.cutoffY; 
  m_eik_params.lambda    = (m_runmode!=run_mode::test)?
    dr->GetValue<double>("lambda",0.3):0.;
  m_eik_params.Delta     = dr->GetValue<double>("Delta",0.4);
  m_eik_params.beta02    = m_ff_params.beta02;
  m_eik_params.bmax      = 2.*m_bmax;
  m_eik_params.accu      = m_accu;
}

void MinBias_Parameters::FillLadderParameters(ATOOLS::Data_Reader * dr) {
  /*
  std::string asf(dr->GetValue<std::string>("As_Form",std::string("smooth")));
  m_as_form = MODEL::asform::smooth;
  if (asf==std::string("constant"))    m_as_form = MODEL::asform::constant;
  else if (asf==std::string("frozen")) m_as_form = MODEL::asform::frozen;
  else if (asf==std::string("smooth")) m_as_form = MODEL::asform::smooth;
  else if (asf==std::string("IR0"))    m_as_form = MODEL::asform::IR0;
  else if (asf==std::string("GDH"))    m_as_form = MODEL::asform::GDH_inspired;
 

  std::string ladderweight = 
    dr->GetValue<std::string>("Ladder_Weight",std::string("Regge"));
  if (ladderweight==std::string("IntervalOnly"))
    m_ladderweight = ladder_weight::IntervalOnly;
  else if (ladderweight==std::string("ReggeDiffusion"))
    m_ladderweight = ladder_weight::ReggeDiffusion;
  else
    m_ladderweight = ladder_weight::Regge;

  std::string ktform(dr->GetValue<std::string>("KT_Form",
					       std::string("smooth")));
  m_ktform = ktform::smooth;
  if (ktform==std::string("cut"))         m_ktform = ktform::cut;
  else if (ktform==std::string("IR0"))    m_ktform = ktform::IR0;
  else if (ktform==std::string("frozen")) m_ktform = ktform::frozen;
  else if (ktform==std::string("smooth")) m_ktform = ktform::smooth;

  std::string ordering =
    dr->GetValue<std::string>("Ordering",std::string("ao_phys"));
  m_ordering = ordering::ao_phys;
  if (ordering==std::string("rap_phys"))      m_ordering = ordering::rap_phys;
  else if (ordering==std::string("ao_phys"))  m_ordering = ordering::ao_phys;
  else if (ordering==std::string("ao_keep"))  m_ordering = ordering::ao_keep;
  else if (ordering==std::string("ao"))       m_ordering = ordering::ao;
  else if (ordering==std::string("keep"))     m_ordering = ordering::keep;
  else if (ordering==std::string("rap_only")) m_ordering = ordering::rap_only;

  std::string resktmin =
    dr->GetValue<std::string>("Resc_KTMin",std::string("off"));
  if (resktmin==std::string("props"))
    m_resc_ktmin = resc_ktmin::props;
  else if (resktmin==std::string("on"))
    m_resc_ktmin = resc_ktmin::on;
  else
    m_resc_ktmin = resc_ktmin::off;

  std::string rescnosing =
    dr->GetValue<std::string>("Resc_NoSinglet",std::string("on"));
  if (rescnosing==std::string("off"))
    m_resc_nosing = resc_nosing::off;
  else
    m_resc_nosing = resc_nosing::on;

  std::string rescoversing =
    dr->GetValue<std::string>("RescOverSinglet",std::string("on"));
  if (rescoversing==std::string("off"))
    m_resc_over_sing = resc_over_sing::off;
  else
    m_resc_over_sing = resc_over_sing::on;


  std::string rescmode =
    dr->GetValue<std::string>("Rescattering",std::string("same"));
  if (rescmode==std::string("off")) {
    m_rescmode = resc_mode::off;
  }
  else {
    m_rescmode  = resc_mode::on;
  }
  m_params["NLaddersFix"] = dr->GetValue<int>("N_Ladders_Fix",-1);
  m_params["KTMin_Mode"]  = dr->GetValue<int>("KTMin_Mode",0);
  m_params["Q02"]         = dr->GetValue<double>("Q_0^2", 3.021183);
  m_params["Q_as2"]       = dr->GetValue<double>("Q_as^2",1.0);
  m_params["Q12"]         = dr->GetValue<double>("Q_1^2",0.0);
  m_params["QN2"]         = dr->GetValue<double>("Q_N^2",0.0);
  m_params["SingletWt"]   = dr->GetValue<double>("Chi_S",0.652533);
  m_params["Ddiff2"]      = dr->GetValue<double>("D_diff^2",0.0);
  m_params["kdiff"]       = dr->GetValue<double>("K_diff",0.0);
  */
}

void MinBias_Parameters::FillShowerLinkParameters(ATOOLS::Data_Reader * dr) {
  /*
  m_params["shower_mode"] = dr->GetValue<int>("Shower_Mode",3);
  m_params["min_kt2"]     = dr->GetValue<double>("Shower_Min_KT2",1.196473);
  m_params["kt2_factor"]  = dr->GetValue<double>("KT2_Factor",3.487842);
  m_params["diff_factor"] = dr->GetValue<double>("Diff_Factor",1.0);
  // rescatterings
  m_params["RescProb"]    = dr->GetValue<double>("RescProb",1.018779);
  m_params["RescProb1"]   = dr->GetValue<double>("RescProb1",0.185517);
  m_params["QRC2"]        = dr->GetValue<double>("Q_RC^2",0.500756);
  m_params["ReconnProb"]  = dr->GetValue<double>("ReconnProb",-15.301277);
  m_params["Misha"]       = dr->GetValue<int>("Misha",0);
  std::string reconnmode =
    dr->GetValue<std::string>("Reconnections",std::string("fix"));
  if (reconnmode==std::string("off")) {
    m_reconnmode = reconn_mode::off;
  }
  else if (reconnmode==std::string("fix")) {
    m_reconnmode = reconn_mode::fix;
  }
  else m_reconnmode = reconn_mode::run;
  */
}


