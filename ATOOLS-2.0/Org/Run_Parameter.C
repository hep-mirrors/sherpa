#include <iostream> 
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Message.H"

using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;

Run_Parameter AORGTOOLS::rpa;

void Run_Parameter::Init(std::string _path,int argc,char* argv[])
{
  path = _path;
  gen.output = 0;
  me.m_string = String_Type::String;
    
  // read in parameters
  Data_Read dr(path+std::string("/Run.dat"));

  gen.output             = dr.GetValue<int>("OUTPUT");
  msg.Init(gen.Output());

  // note: "particle_init" must have been called before calling these functions!
  gen.beam1              = dr.GetValue<Flavour>("BEAM1");
  gen.beam2              = dr.GetValue<Flavour>("BEAM2");

  gen.nevents            = dr.GetValue<long>("EVENTS");
  gen.cms_energy         = dr.GetValue<double>("CMSENERGY");
  gen.accu               = dr.GetValue<double>("Num. Accuracy");
  gen.runtime            = dr.GetValue<std::string>("Runtime"); // Time
  gen.massive            = dr.GetValue<Switch::code>("MASSES");

  me.model_file          = dr.GetValue<std::string>("MODELFILE");
  me.model_data_file     = dr.GetValue<std::string>("MODEL_DATA_FILE");
  me.the_model           = dr.GetValue<Model_Type::code>("MODEL");
  me.model_mass          = dr.GetValue<Switch::code>("MODELMASS");
  me.m_string            = dr.GetValue<String_Type::code>("CALCULATION");
  if (me.m_string==NotDefined<String_Type::code>())
    me.m_string            = dr.GetValue<String_Type::code>("STRING");

  me.cut_scheme          = dr.GetValue<int>("CUTSCHEME");
  if (me.cut_scheme==NotDefined<int>()) 
    me.cut_scheme = 0;

  me.coulomb_corr        = dr.GetValue<Switch::code>("COULOMB");
  me.kfactorscheme       = dr.GetValue<int>("KFACTORSCHEME");
  me.scalescheme         = dr.GetValue<int>("SCALESCHEME");
  
  // read in some consts
  dr.ReadIn(path+std::string("/")+me.model_file);

  // read in comand line
  int iarg=0;
  while ((argc-iarg)>=3) {
    if (argv[iarg][0]=='-' && argv[iarg][1]=='p') {
      dr.SetValue(std::string(argv[iarg+1]),std::string(argv[iarg+2]));
      iarg+=3;
    }
    else ++iarg;
  }

  consts.aqed_mz         = 1./dr.GetValue<double>("alpha_QED(MZ)");
  consts.sin2_tw         = dr.GetValue<double>("SIN2_TW");
  if (consts.sin2_tw==NotDefined<double>())
    consts.sin2_tw         = dr.GetValue<double>("sinTW^2");
  consts.running_masses  = dr.GetValue<Switch::code>("RUNMASS");
  consts.running_widths  = dr.GetValue<int>("RUNWIDTH");
  consts.running_aqed    = dr.GetValue<Switch::code>("RUNAQED");
  consts.running_alphas  = dr.GetValue<Switch::code>("RUN_ALPHA_S");

  consts.alpha_s_fix     = dr.GetValue<double>("ALPHA_S_FIX");
  if (consts.alpha_s_fix==NotDefined<double>())
    consts.alpha_s_fix     = dr.GetValue<double>("alpha_S_fixed");
  if (consts.alpha_s_fix==NotDefined<double>()) {
    msg.Out()<<" WARNING using alpha_S(MZ) for alpha_S_fixed "<<std::endl;
    consts.alpha_s_fix     = dr.GetValue<double>("alpha_S(MZ)"); 
  }

  pshower.final_q02      = dr.GetValue<double>("FINAL_Q02");
  pshower.lambda_qcd     = dr.GetValue<double>("LAMBDA_QCD"); //(4) or (5) ?!!!
  pshower.initial_q02    = - dr.GetValue<double>("INITIAL_Q02"); // note sign change

  dshower.k_perp_scheme  = dr.GetValue<int>("K_PERP_SCHEME");
  dshower.k_ord_scheme   = dr.GetValue<int>("K_ORD_SCHEME");
  dshower.kt_min         = dr.GetValue<double>("K_T_MIN");

  test.analysis          = dr.GetValue<int>("ANALYSIS");

  test.facycut = dr.GetValue<double>("TEST_FAC_YCUT");
  if (test.facycut==NotDefined<double>()) test.facycut = 0.125;
  test.facnlly = dr.GetValue<double>("TEST_FAC_NLLQ");
  if (test.facnlly==NotDefined<double>()) test.facnlly = 1.;


  test.faca = dr.GetValue<double>("TEST_FAC_A");
  if (test.faca==NotDefined<double>()) test.faca = 1.;
  test.facb = dr.GetValue<double>("TEST_FAC_B");
  if (test.facb==NotDefined<double>()) test.facb = 1.;
  test.facc = dr.GetValue<double>("TEST_FAC_C");
  if (test.facc==NotDefined<double>()) test.facc = 1.; 

  test.flaga = dr.GetValue<int>("TEST_FLAG_A");
  if (test.flaga==NotDefined<int>()) test.flaga = 0;
  test.flagb = dr.GetValue<int>("TEST_FLAG_B");
  if (test.flagb==NotDefined<int>()) test.flagb = 0;
  test.flagc = dr.GetValue<int>("TEST_FLAG_C");
  if (test.flagc==NotDefined<int>()) test.flagc = 0; 

  dr.ReadIn(path+std::string("/Integration.dat"));

  integ.ycut             = dr.GetValue<double>("YCUT");
  integ.accu             = dr.GetValue<double>("ERROR");
  integ.integrator       = dr.GetValue<int>("INTEGRATOR");
  integ.jetfinder        = dr.GetValue<int>("JETFINDER");

  gen.rpa_id = dr.GenerateKey();
  //  dr.WriteOut("current_flags.log");
  //  dr.WriteOut(std::string("save/rpa/")+gen.rpa_id+std::string(".dat"));

  if(gen.Masses()==Switch::Off) APHYTOOLS::SetMassless();
}

Run_Parameter::~Run_Parameter() { }


double Run_Parameter::Consts::Mass(Flavour flav,double t)
{
  if (flav.IsBoson()) {
    if (IsWidthRunning()==1) {
      double gamm = flav.Width()/flav.Mass();
      return flav.Mass()/sqrt(1.-gamm*gamm);
    } 
    return flav.Mass();
  }
  if (rpa.gen.Masses()==Switch::Off) return 0.;
  return flav.Mass();
}  

double Run_Parameter::Consts::Width(Flavour flav,double t)
{
  if (flav.IsBoson()) {
    if (IsWidthRunning() == 1) {
      double gamm = flav.Width()/flav.Mass();
      return flav.Width()/sqrt(1.-gamm*gamm);
    } 
    if (IsWidthRunning() == 2) return flav.Width()*sqrt(t)/flav.Mass();
    // *AS*  in Zfitter:     if (Run_Width() == 2) return flav.Width()*t/sqr(flav.Mass());
    return flav.Width();
  }
  //later on running width for fermions like top
  return flav.Width();
}

double Run_Parameter::Consts::RunningMass(Flavour flav,double t)
{
  if (flav.Mass()==0.) return 0.;

  if (flav.IsBoson()) {
    if (IsWidthRunning() == 1) {
      double gamm = flav.Width()/flav.Mass();
      return flav.Mass()/sqrt(1.-gamm*gamm);
    } 
    return flav.Mass();
  
  }
  // turn off running by uncommenting:
  //  return flav.Mass();


  // test 0.3 GeV limit for masses
  //*AS*    if (sqr(flav.Mass()) < 0.3) return 0.;
  if (t<1.) t=1.;
  if ((t<sqr(flav.Mass())) || (!flav.IsQuark())) return flav.Mass(); 
}  

