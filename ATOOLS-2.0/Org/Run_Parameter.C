#include <iostream> 
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Message.H"

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
  gen.p_model     = NULL;
} 

void Run_Parameter::Init(std::string path,std::string file,int argc,char* argv[])
{
  m_path        = path;
  Data_Read dr(m_path+file);
  gen.m_output             = dr.GetValue<int>("OUTPUT",0);
  gen.m_analysis           = dr.GetValue<int>("ANALYSIS",0);
  gen.m_nevents            = dr.GetValue<long>("EVENTS",100);
  gen.m_accu               = dr.GetValue<double>("Num. Accuracy",1.e-10);
  //gen.m_runtime            = dr.GetValue<std::string>("Runtime"); // Time

  msg.Init(gen.Output());
  gen.m_rpa_id = dr.GenerateKey();
}

Run_Parameter::~Run_Parameter() { }


