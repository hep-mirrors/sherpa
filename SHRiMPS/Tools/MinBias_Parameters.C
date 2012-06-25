#include "SHRiMPS/Tools/MinBias_Parameters.H" 
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;




MinBias_Parameters SHRIMPS::MBpars;

MinBias_Parameters::MinBias_Parameters() {}

MinBias_Parameters::~MinBias_Parameters() {}

void MinBias_Parameters::Init(ATOOLS::Data_Reader * dr) {
  // general properties: c.m. energy in rapidity and its cutoff, 
  // impact parameters
  m_params["originalY"]   = log(ATOOLS::rpa->gen.Ecms()/
			       ATOOLS::Flavour(kf_p_plus).HadMass());
  m_params["deltaY"]      = dr->GetValue<double>("deltaY",1.62);
  m_params["bmin"]        = dr->GetValue<double>("bmin",0.);
  m_params["bmax"]        = dr->GetValue<double>("bmax",20.);
  m_params["accu"]        = dr->GetValue<double>("accu",5.e-4);
  // form factors
  m_params["NGWstates"]   = dr->GetValue<int>("GW_States",2);
  m_params["FFpref"]      = 1./sqrt(m_params["NGWstates"]);
  m_params["Lambda2"]     = dr->GetValue<double>("Lambda2",0.89);
  m_params["beta02(mb)"]  = dr->GetValue<double>("beta_0^2",25.0);
  m_params["beta0"]       = sqrt(1.e9*m_params["beta02(mb)"]/
				 ATOOLS::rpa->Picobarn());
  m_params["kappa"]       = dr->GetValue<double>("kappa",0.36);
  m_params["xi"]          = dr->GetValue<double>("xi",0.01);
  // parameters of the eikonal
  m_params["lambda"]      = dr->GetValue<double>("lambda",0.56);
  m_params["Delta"]       = dr->GetValue<double>("Delta",0.28);
  // ladder generation
  m_params["NLaddersFix"] = dr->GetValue<int>("N_Ladders_Fix",-1);
  m_params["KTMin_Mode"]  = dr->GetValue<int>("KTMin_Mode",0);
  m_params["Q02"]         = dr->GetValue<double>("Q_0^2",1.);
  m_params["Q_as2"]       = dr->GetValue<double>("Q_as^2",m_params["Q02"]);
  m_params["Q12"]         = dr->GetValue<double>("Q_1^2",0.);
  m_params["QN2"]         = dr->GetValue<double>("Q_N^2",0.);
  m_params["SingletWt"]   = dr->GetValue<double>("Chi_S",1.);
  m_params["Ddiff2"]      = dr->GetValue<double>("D_diff^2",0.);
  m_params["kdiff"]       = dr->GetValue<double>("K_diff",0.);
  // showering off soft stuff
  m_params["shower_mode"] = dr->GetValue<int>("Shower_Mode",3);
  m_params["min_kt2"]     = dr->GetValue<double>("Shower_Min_KT2",16.);
  m_params["kt2_factor"]  = dr->GetValue<double>("KT2_Factor",4.);
  m_params["diff_factor"] = dr->GetValue<double>("Diff_Factor",1.);
  // rescatterings
  m_params["RescProb"]    = dr->GetValue<double>("RescProb",1.);
  m_params["RescProb1"]   = dr->GetValue<double>("RescProb1",0.);
  m_params["SpatProb"]    = dr->GetValue<double>("SpatProb",0.);
  m_params["SpatWidth"]   = dr->GetValue<double>("SpatWidth",
						 m_params["Lambda2"]);
  m_params["QRC2"]        = dr->GetValue<double>("Q_RC^2",4.*m_params["Q02"]);
  m_params["ReconnProb"]  = dr->GetValue<double>("ReconnProb",1.);

  std::string ffform = 
    dr->GetValue<std::string>("FF_Form",std::string("dipole"));
  if (ffform==std::string("dipole"))
    m_ffform = ff_form::dipole;
  else
    m_ffform = ff_form::Gauss;

  std::string absorption = 
    dr->GetValue<std::string>("Absorption",std::string("factorial"));
  if (absorption==std::string("exponential"))
    m_absorp = absorption::exponential;
  else
    m_absorp = absorption::factorial;

  std::string asf(dr->GetValue<std::string>("As_Form",std::string("IR0")));
  m_as_form = MODEL::asform::smooth;
  if (asf==std::string("constant"))    m_as_form = MODEL::asform::constant;
  else if (asf==std::string("frozen")) m_as_form = MODEL::asform::frozen;
  else if (asf==std::string("smooth")) m_as_form = MODEL::asform::smooth;
  else if (asf==std::string("IR0"))    m_as_form = MODEL::asform::IR0;
  else if (asf==std::string("GDH"))    m_as_form = MODEL::asform::GDH_inspired;
 
  std::string runmode =
    dr->GetValue<std::string>("Shrimps_Mode",std::string("Inelastic"));
  if (runmode==std::string("Xsecs")) 
    m_runmode = run_mode::xsecs_only;
  else if (runmode==std::string("Elastic")) 
    m_runmode = run_mode::elastic_events;
  else if (runmode==std::string("Single-diffractive")) 
    m_runmode = run_mode::single_diffractive_events;
  else if (runmode==std::string("Double-diffractive")) 
    m_runmode = run_mode::double_diffractive_events;
  else if (runmode==std::string("Quasi-elastic")) 
    m_runmode = run_mode::quasi_elastic_events;
  else if (runmode==std::string("Inelastic")) 
    m_runmode = run_mode::inelastic_events;
  else if (runmode==std::string("All")) 
    m_runmode = run_mode::all_min_bias;
  else if (runmode==std::string("Underlying")) 
    m_runmode = run_mode::underlying_event;

  std::string weightmode =      
    dr->GetValue<std::string>("MB_Weight_Mode",std::string("Unweighted"));
  if (weightmode==std::string("Unweighted")) 
    m_weightmode = weight_mode::unweighted;
  else if (weightmode==std::string("Weighted")) 
    m_weightmode = weight_mode::weighted;

  std::string ladderweight = 
    dr->GetValue<std::string>("Ladder_Weight",std::string("Regge"));
  if (ladderweight==std::string("IntervalOnly"))
    m_ladderweight = ladder_weight::IntervalOnly;
  else if (ladderweight==std::string("ReggeDiffusion"))
    m_ladderweight = ladder_weight::ReggeDiffusion;
  else
    m_ladderweight = ladder_weight::Regge;

  std::string ktform(dr->GetValue<std::string>("KT_Form",std::string("IR0")));
  m_ktform = ktform::smooth;
  if (ktform==std::string("cut"))         m_ktform = ktform::cut;
  else if (ktform==std::string("IR0"))    m_ktform = ktform::IR0;
  else if (ktform==std::string("frozen")) m_ktform = ktform::frozen;
  else if (ktform==std::string("smooth")) m_ktform = ktform::smooth;

  std::string ordering =
    dr->GetValue<std::string>("Ordering",std::string("ao_phys"));
  if (ordering==std::string("rap_phys"))
    m_ordering = ordering::rap_phys;
  else if (ordering==std::string("ao_keep"))
    m_ordering = ordering::ao_keep;
  else if (ordering==std::string("ao"))
    m_ordering = ordering::ao;
  else if (ordering==std::string("keep"))
    m_ordering = ordering::keep;
  else if (ordering==std::string("rap_only"))
    m_ordering = ordering::rap_only;
  else
    m_ordering = ordering::rap_only;

  std::string resktmin =
    dr->GetValue<std::string>("Resc_KTMin",std::string("on"));
  if (resktmin==std::string("props"))
    m_resc_ktmin = resc_ktmin::props;
  else if (resktmin==std::string("on"))
    m_resc_ktmin = resc_ktmin::on;
  else
    m_resc_ktmin = resc_ktmin::off;

  std::string rescnosing =
    dr->GetValue<std::string>("Resc_NoSinglet",std::string("off"));
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
  else if (rescmode==std::string("space")) {
    m_rescmode  = resc_mode::on;
    m_spacemode = spat_mode::expo;
  }
  else {
    m_rescmode  = resc_mode::on;
    m_spacemode = spat_mode::off;
  }

  std::string reconnmode =
    dr->GetValue<std::string>("Reconnections",std::string("fix"));
  if (reconnmode==std::string("off")) {
    m_reconnmode = reconn_mode::off;
  }
  else if (reconnmode==std::string("fix")) {
    m_reconnmode = reconn_mode::fix;
  }
  else m_reconnmode = reconn_mode::run;
}

double MinBias_Parameters::operator()(std::string keyword) {
  if (m_params.find(keyword)==m_params.end()) {
    msg_Error()<<"Error in "<<METHOD<<"("<<keyword<<") : "<<std::endl
	       <<"   Keyword not found, will exit the run."<<std::endl;
    abort();
  }
  return m_params[keyword];
}

std::ostream & SHRIMPS::operator<<(std::ostream & s, 
				 const run_mode::code & runmode) {
  switch (runmode) {
  case run_mode::xsecs_only:
    s<<" XSecs only";
    break;
  case run_mode::elastic_events:
    s<<" Elastic events";
    break;
  case run_mode::single_diffractive_events:
    s<<" Single-diffractive events";
    break;
  case run_mode::quasi_elastic_events:
    s<<" Quasi-elastic events";
    break;
  case run_mode::inelastic_events:
    s<<" Inelastic events";
    break;
  case run_mode::all_min_bias:
    s<<" Minimum Bias";
    break;
  case run_mode::underlying_event:
    s<<" Underlying event";
    break;
  default:
    s<<" Unknown mode";
    break;
  }
  return s;
}

std::ostream & SHRIMPS::operator<<(std::ostream & s, 
				 const weight_mode::code & weightmode) {
  switch (weightmode) {
  case weight_mode::weighted:
    s<<" Weighted events";
    break;
  case weight_mode::unweighted:
    s<<" Unweighted events";
    break;
  default:
    s<<" Unknown mode";
    break;
  }
  return s;
}


std::ostream & SHRIMPS::operator<<(std::ostream & s, const ff_form::code & fff) {
  if (fff==ff_form::Gauss)       s<<"Gaussian";
  else if (fff==ff_form::dipole) s<<"Dipole";
  else                           s<<"unknown";
  return s;
}

std::ostream & SHRIMPS::operator<<(std::ostream & s, 
				 const absorption::code & a) {
  if (a==absorption::exponential)    s<<"exponential";
  else if (a==absorption::factorial) s<<"factorial";
  else                               s<<"unknown";
  return s;
}

std::ostream & operator<<(std::ostream & s, 
			  const ladder_weight::code & wt) {
  if (wt==ladder_weight::IntervalOnly)        s<<"y interval only ";
  else if (wt==ladder_weight::ReggeDiffusion) s<<"Regge factor + diffusion ";
  else if (wt==ladder_weight::Regge)          s<<"Regge factor ";
  else                                    s<<"unknown";
  return s;
}
