#include <iostream> 
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Message.H"

#include "Running_AlphaS.H"

using namespace AMATOOLS;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace std;

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
      dr.SetValue(string(argv[iarg+1]),string(argv[iarg+2]));
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
    msg.Out()<<" WARNING using alpha_S(MZ) for alpha_S_fixed "<<endl;
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
  //  dr.WriteOut(string("save/rpa/")+gen.rpa_id+string(".dat"));
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

  /*  
//======================================================================
//
//   Mass running:   MSBAR as in Dittmaier
//
//                                    alphaS(t)
// mass_b(t) = mass_b(mass_b) * ( 1 - -------- ( ln (t/mass_b^2) + 4/3))
//                                      pi
//

  double mass0   = flav.Mass();
  double as_t    = mo->Aqcd(t);

  //  for (double s=6. ; s<120; s+=5)
  //    cout<<" s="<<s <<" as="<<mo->Aqcd(sqr(s))<<endl;

  double mass = mass0 * ( 1. - as_t / M_PI * ( log( t/sqr(mass0)) + 4./3.));
  cout<<" in Run_Parameter::RunningMass() "<<endl;
  cout<<" t="<<t<<" mass0= "<<mass0<<" => "<<mass<<endl;

  return mass;
  */

  //=====================================================================
  //  if (Run_Mass()==Switch::On) {

    // running for quarks only
    // alpha_S(M_Z^2)
    double asMZ=(*as)(sqr(Flavour(kf::Z).PSMass()));
    //    cout<<" in Run_Parameter::Mass() "<<endl;
    //    cout<<" alpha_s(Mz^2)="<<asMZ<<endl;
    // alpha_S(t)
    double as2=(*as)(t);
    //    cout<<" alpha_s("<<t<<")="<<as2<<endl;

    double factor = 1.;
    // expo = 6 * CF/(11 - 2/3 n_f)     /2 !!!!
    double mas_Z  = Flavour(kf::Z).PSMass();
    double mas_s2 = Flavour(kf::s).PSMass()*Flavour(kf::s).Mass();
    double mas_c2 = Flavour(kf::c).PSMass()*Flavour(kf::c).Mass();
    double mas_b2 = Flavour(kf::b).PSMass()*Flavour(kf::b).Mass();
    double mas_t2 = Flavour(kf::t).PSMass()*Flavour(kf::t).Mass();
    
    double ast = (*as)(mas_t2);
    double asb = (*as)(mas_b2);
    double asc = (*as)(mas_c2);
    double ass = (*as)(mas_s2);

    if (t>mas_t2) {
      factor*= pow(ast/asb,4./(11.-10./3.));
      factor*= pow(as2/ast,4./(11.-12./3.));
    }
    else if (t>mas_b2) factor*= pow(as2/asb,4./(11.-10./3.));
    if (((flav==Flavour(kf::c)) || flav==Flavour(kf::c).Bar()) && (t>mas_b2)) 
      factor*= pow(asb/asc,4./(11.-8./3.));
    if (((flav==Flavour(kf::c)) || flav==Flavour(kf::c).Bar())&& (t<=mas_b2)) 
      factor*= pow(as2/asc,4./(11.-8./3.));  
    // running quark masses below 0.3 GeV are assumed to be massless!!!
    if (sqr(flav.Mass()*factor) < 0.3) factor = 0.;

    double as_M =(*as)(sqr(flav.Mass()));
    double factor_b=1.- 4. * as_M / (3. * M_PI);

    factor*=factor_b; //1.- 4. * as2 / (3. * M_PI);
    return flav.Mass()*factor;
    // }
  return flav.Mass();
}  

