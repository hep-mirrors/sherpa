#include <iostream> 
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Message.H"
#include "Random.H"

using namespace ATOOLS;

Run_Parameter ATOOLS::rpa;

Run_Parameter::Run_Parameter() 
{
  gen.m_output    = gen.m_analysis   = 0;
  gen.m_nevents   = 0;
  gen.m_cutscheme = 0;
  gen.m_ecms      = gen.m_ycut       = gen.m_accu        = 0.;
  gen.m_beam1     = gen.m_beam2      = Flavour(kf::none);
  gen.m_rpa_id    = gen.m_flavour_id = std::string("");
#ifndef USING__ATOOLS_only
  gen.p_model     = NULL;
#endif
} 

void Run_Parameter::Init(std::string path,std::string file,int argc,char* argv[])
{
  gen.m_timer.Start();
  std::string gccversion;
  system("gcc -dumpversion > sherpa_gcc_test");
  std::ifstream *gcctest = new std::ifstream("sherpa_gcc_test");
  if (*gcctest) (*gcctest)>>gccversion;
  delete gcctest;
  system("if test -f sherpa_gcc_test; then rm sherpa_gcc_test; fi");
  if (gccversion==std::string("2.96")) {
    ATOOLS::msg.Error()<<"Run_Parameter::Init(..): "<<om::red
		       <<"Sherpa must not be run with gcc version 2.96 !"<<om::reset<<std::endl;
    exit(1);
  }
  m_path        = path;
  Data_Read dr(m_path+file);
  gen.m_output             = dr.GetValue<int>("OUTPUT",0);
  std::string logfile=dr.GetValue<std::string>("LOG_FILE",std::string(""));
  msg.Init(gen.Output(),logfile);
  gen.m_analysis           = dr.GetValue<int>("ANALYSIS",0);
  gen.m_nevents            = dr.GetValue<long>("EVENTS",100);
  // read only if defined (no error message if not defined)
  gen.m_seed               = dr.GetValue<long>("RANDOM_SEED");
  if (gen.m_seed==NotDefined<long>()) gen.m_seed=1234;
  gen.m_timeout            = dr.GetValue<double>("TIMEOUT");
  if (gen.m_timeout<0.) gen.m_timeout=0.;
  rpa.gen.m_timer.Start();
  gen.m_batchmode          = dr.GetValue<int>("BATCH_MODE");
  if (gen.m_batchmode==NotDefined<int>()) gen.m_batchmode=1;
  double ycut=dr.GetValue<double>("YCUT");
  if (ycut!=NotDefined<double>()) gen.m_ycut=ycut;
  gen.m_accu               = dr.GetValue<double>("Num. Accuracy",1.e-10);
  //gen.m_runtime            = dr.GetValue<std::string>("Runtime"); // Time
  Switch::code color=dr.GetValue<Switch::code>("PRETTY_PRINT");
  if (color==Switch::On) msg.SetModifiable(true);
  gen.m_rpa_id = dr.GenerateKey();
  if (gen.m_seed!=1234) ran.SetSeed(gen.m_seed);
}

Run_Parameter::~Run_Parameter() 
{ 
  if (msg.Level()>=1) gen.m_timer.PrintTime();
}

bool Run_Parameter::Gen::CheckTime(const double limit)
{ 
  if (limit==0.) {
    if (m_timeout>0.) return m_timer.UserTime()<m_timeout;
  }
  else {
    return m_timer.UserTime()<limit;
  }
  return false;
}


