#include <iostream> 
#include "Run_Parameter.H"
#include "MathTools.H"
#include "Message.H"
#include "Exception.H"
#include "Random.H"
#include "Data_Reader.H"
#include "Data_Collector.H"
#include <stdlib.h>

using namespace ATOOLS;

bool ATOOLS::Run_Parameter::s_initialized=false;
std::map<std::string,std::string> ATOOLS::Run_Parameter::s_variables;
Run_Parameter ATOOLS::rpa;

Run_Parameter::Run_Parameter() 
{
  AnalyseEnvironment();
  gen.m_output    = gen.m_analysis   = 0;
  gen.m_nevents   = 0;
  gen.m_cutscheme = 0;
  gen.m_ecms      = gen.m_ycut       = gen.m_accu        = 0.;
  gen.m_delta_r   = 1.;
  gen.m_beam1     = gen.m_beam2      = Flavour(kf::none);
  gen.m_rpa_id    = gen.m_flavour_id = std::string("");
#ifndef USING__ATOOLS_only
  gen.p_model     = NULL;
#endif
  gen.m_ndicedevents = 0;
} 

void Run_Parameter::AnalyseEnvironment() 
{
  if (s_initialized) return;
#ifdef __GNUC__
#if __GNUC__ == 2 && __GNUC_MINOR__ == 96
#error Sherpa was not designed for gcc 2.96
#endif
#endif
  std::string gccversion;
  system("gcc -dumpversion > sherpa_gcc_test");
  std::ifstream *test = new std::ifstream("sherpa_gcc_test");
  if (*test) (*test)>>gccversion;
  delete test;
  system("if test -f sherpa_gcc_test; then rm sherpa_gcc_test; fi");
  if (gccversion.find("2.96")!=std::string::npos) {
    THROW(fatal_error,"Sherpa must not be run on gcc version 2.96 !");
  }
  char *var=NULL;
  s_variables["PATH"]=std::string(((var=getenv("PATH"))==NULL?"":var));
  s_variables["SHERPASYS"]=std::string(((var=getenv("SHERPASYS"))==NULL?"":var));
  s_variables["SHERPA_PDF_PATH"]=std::string(((var=getenv("SHERPA_PDF_PATH"))==NULL?"":var));
  s_variables["SHERPA_CPP_PATH"]=std::string(((var=getenv("SHERPA_CPP_PATH"))==NULL?"":var));
  s_variables["SHERPA_LIB_PATH"]=std::string(((var=getenv("SHERPA_LIB_PATH"))==NULL?"":var));
  s_variables["LD_LIBRARY_PATH"]=std::string(((var=getenv("LD_LIBRARY_PATH"))==NULL?"":var));
  if (system("test -f Sherpa")) {
    std::string paths=s_variables["PATH"];
    do {
      size_t pos=ATOOLS::Min(paths.length(),paths.find(":"));
      std::string cur=paths.substr(0,pos);
      if (!system((std::string("test -f ")+cur+std::string("/Sherpa")).c_str())) {
	s_variables["SHERPA_BIN_PATH"]=cur;
	break;
      }
      if (pos<paths.length()) paths=paths.substr(pos+1);
      else paths="";
    } while (paths.length()>0);
  }
  s_variables["SHERPA_RUN_PATH"]=getenv("PWD");
  s_initialized=true;
}

void Run_Parameter::Init(std::string path,std::string file,int argc,char* argv[])
{
  gen.m_timer.Start();
  system("finger `whoami` > sherpa_user_test");
  Data_Reader *reader = new Data_Reader();
  reader->SetInputFile("sherpa_user_test");
  reader->SetVectorType(reader->VHorizontal);
  std::vector<std::string> help;
  if (!reader->VectorFromFile(help,"Name:")) { 
    gen.m_username=std::string("<unknown user>");
  }
  else {
    for (std::vector<std::string>::iterator nit=help.begin();nit!=help.end();++nit) {
      gen.m_username+=*nit+std::string(" ");
    }
  }
  delete reader;
  msg.Out()<<"Welcome to Sherpa, "<<gen.m_username
	   <<". Initialization of framework underway."<<std::endl;
  system("if test -f sherpa_user_test; then rm sherpa_user_test; fi");
  m_path = path;
  Data_Read dr(m_path+file);
  gen.m_output = dr.GetValue<int>("OUTPUT",0);
  std::string logfile=dr.GetValue<std::string>("LOG_FILE",std::string(""));
  msg.Init(gen.m_output,logfile);
  // make path nice
  if (path.length()>0) {
    if (path[0]!='/') path=std::string(getenv("PWD"))+std::string("/")+path;
    while (path[path.length()-1]=='/' || path[path.length()-1]=='.') 
      path=path.substr(0,path.length()-1);
  }
  // search for sherpa executeable
  if (argc>0) {
    std::string command=argv[0];
    if (!system((std::string("test -f ")+command).c_str())) {
      if (command[0]!='/') command=std::string(getenv("PWD"))+std::string("/Sherpa");
      s_variables["SHERPA_BIN_PATH"]=command.substr(0,command.length()-7);
    }
    else if (!system("test -f $PWD/Sherpa")) {
      s_variables["SHERPA_BIN_PATH"]=std::string(getenv("PWD"));
    }
  }
  // set pdf path
  std::string pdfpath=dr.GetValue<std::string>("SHERPA_PDF_PATH",std::string(""));
  if (pdfpath.length()>0 && pdfpath[0]=='/') s_variables["SHERPA_PDF_PATH"]=pdfpath;
  else if (s_variables["SHERPA_PDF_PATH"].length()==0) 
    s_variables["SHERPA_PDF_PATH"]=s_variables["SHERPA_BIN_PATH"];
  // set cpp path
  std::string cpppath=dr.GetValue<std::string>("SHERPA_CPP_PATH",std::string(""));
  if (cpppath.length()>0 && cpppath[0]=='/') s_variables["SHERPA_CPP_PATH"]=cpppath;
  else if (path!=getenv("PWD")) s_variables["SHERPA_CPP_PATH"]=path;
  else if (s_variables["SHERPA_CPP_PATH"].length()==0) 
    s_variables["SHERPA_CPP_PATH"]=s_variables["SHERPA_RUN_PATH"];
  // set lib path
  std::string libpath=dr.GetValue<std::string>("SHERPA_LIB_PATH",std::string(""));
  if (libpath.length()>0 && libpath[0]=='/') s_variables["SHERPA_LIB_PATH"]=libpath;
  else if (s_variables["SHERPA_LIB_PATH"].length()==0) 
    s_variables["SHERPA_LIB_PATH"]=s_variables["SHERPA_CPP_PATH"]
      +std::string("/Process/lib");
  msg_Tracking()<<"Run_Parameter::Init(..): Paths are {\n"
		<<"   SHERPA_BIN_PATH = "<<s_variables["SHERPA_BIN_PATH"]<<"\n"
		<<"   SHERPA_PDF_PATH = "<<s_variables["SHERPA_PDF_PATH"]<<"\n"
		<<"   SHERPA_CPP_PATH = "<<s_variables["SHERPA_CPP_PATH"]<<"\n"
		<<"   SHERPA_LIB_PATH = "<<s_variables["SHERPA_LIB_PATH"]<<"\n"
		<<"}"<<std::endl;
  s_variables["CURRENT_SHERPASYS"]=s_variables["SHERPA_BIN_PATH"]+std::string("/../..");

#ifndef __sgi
  setenv("LD_LIBRARY_PATH",(s_variables["LD_LIBRARY_PATH"]+std::string(":")+
			    s_variables["SHERPA_LIB_PATH"]).c_str(),1);
#endif
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
  int stacktrace           = dr.GetValue<int>("STACK_TRACE");
  if (stacktrace==NotDefined<int>()) stacktrace=0;
  Exception_Handler::SetStackTrace(stacktrace);
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

void  Run_Parameter::Gen::SetEcms(double _ecms)     { 
  m_ecms    = _ecms;
  Data_Collector::AddData("ECMS",new Blob_Data<double>(m_ecms));
}
void  Run_Parameter::Gen::SetBeam1(const Flavour b) { 
  m_beam1  = b;   
  Data_Collector::AddData("BEAM1",new Blob_Data<std::string>(m_beam1.TexName()));
}
void  Run_Parameter::Gen::SetBeam2(const Flavour b) { 
  m_beam2  = b;   
  Data_Collector::AddData("BEAM2",new Blob_Data<std::string>(m_beam2.TexName()));
}

