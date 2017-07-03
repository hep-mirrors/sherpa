#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>
#include <sys/stat.h>

#include "GoSam_Interface.H"

using namespace GoSam;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


extern "C" {

  // initialisation
  void OLP_Start(const char * filename, int* success);
  void gosam_start()
  {
    int success(0);
    OLP_Start("none",&success);
  }

  // finalisation
  void gosam_finish() {}

  // information
  void OLP_Info(char * olp_name,char * olp_version,char * message);
  void gosam_info(char* str)
  { OLP_Info(str,str,str); }

  // parameter setting
  void OLP_SetParameter(const char* name, const double* re, const double* im,
                        int* ierr);
  void gosam_setparameter_complex(const char* key, const double* reval,
                                  const double* imval, int* ierr)
  { OLP_SetParameter(key,reval,imval,ierr); }
  void gosam_setparameter_double(const char* key, const double* val, int* ierr)
  { const double im(0.); OLP_SetParameter(key,val,&im,ierr); }
  void gosam_setparameter_int(const char* key, const int* val, int* ierr)
  { const double re(*val), im(0.); OLP_SetParameter(key,&re,&im,ierr); }
  void gosam_setparameter_string(const char* key, const char* val, int* ierr)
  { msg_Out()<<METHOD<<"(): Nothing done. "<<key<<" "<<val<<std::endl; *ierr=2; }

  // process handling
  void OLP_GetProcessnumber(const char* process, int * nr);
  int  gosam_register_process(const char* process, int * newlibs)
  {
    int nr(-1);
    OLP_GetProcessnumber(process,&nr);
    if (nr<0) THROW(fatal_error,"Process not found.");
    return nr;
  }

  // calculation
  void OLP_EvalSubProcess(int* id, double* pp, double* mu, double* as, double* rval);
  void OLP_EvalSubProcess2(int* id, double* pp, double* mu, double* rval, double* acc);
  void OLP_EvalSubProcess_EW(int* id, double* pp, double* mu, double* rval, double* acc);
  void gosam_evaluate_loop(int* id, double* pp, double* mu, double* rval, double* acc)
  {
    OLP_EvalSubProcess_EW(id,pp,mu,rval,acc);
  }
}


std::map<int,std::string> GoSam_Interface::s_procmap;

std::string GoSam_Interface::s_gosamprefix     = std::string("");
bool        GoSam_Interface::s_ignore_model    = false;
bool        GoSam_Interface::s_exit_on_error   = true;
size_t      GoSam_Interface::s_vmode           = 0;
bool        GoSam_Interface::s_newlibs         = false;

GoSam_Interface::GoSam_Interface() :
  ME_Generator_Base("GoSam")
{}

GoSam_Interface::~GoSam_Interface()
{
  gosam_finish();
}

bool GoSam_Interface::Initialize(const string &path,const string &file,
                                     MODEL::Model_Base *const model,
                                     BEAM::Beam_Spectra_Handler *const beam,
                                     PDF::ISR_Handler *const isr)
{
  // find GoSam installation prefix with several override options
  struct stat st;
  Data_Reader reader(" ",";","#","=");
  s_ignore_model = reader.GetValue<int>("GOSAM_IGNORE_MODEL",0);
  if (s_ignore_model) msg_Info()<<METHOD<<"(): GoSam will use the "
                                <<"Standard Model, even if you set a "
                                <<"different model, without warning."
                                <<std::endl;
  s_exit_on_error = reader.GetValue<int>("GOSAM_EXIT_ON_ERROR",1);
  s_gosamprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/GoSam";
  if(stat(s_gosamprefix.c_str(),&st) != 0) s_gosamprefix = GOSAM_PREFIX;
  s_gosamprefix = reader.GetValue<string>("GOSAM_PREFIX", s_gosamprefix);
  msg_Info()<<"Initialising GoSam generator from "<<s_gosamprefix<<endl;
  s_vmode = reader.GetValue<int>("GOSAM_VMODE", 0);
  msg_Tracking()<<METHOD<<"(): Set V-mode to "<<s_vmode<<endl;

  // load library dynamically
  s_loader->AddPath(s_gosamprefix+"/lib");
  if (!s_loader->LoadLibrary("gosam")) THROW(fatal_error, "Failed to load libgosam.");

  // set OL verbosity
  std::string gosam_verbosity = reader.GetValue<std::string>("GOSAM_VERBOSITY","0");
  SetParameter("Verbosity",gosam_verbosity);

  // tell OL about the current model and check whether accepted
  if (!s_ignore_model) SetParameter("Model", MODEL::s_model->Name());

  // we give parameters to GoSam as alpha(MZ) and masses
  SetParameter("EWScheme","alphaMZ");
  // ew-renorm-scheme to Gmu by default
  switch (ToType<int>(rpa->gen.Variable("EW_REN_SCHEME"))) {
  case 1:
    SetParameter("EWRenormalisationScheme","alpha0");
    break;
  case 2:
    SetParameter("EWRenormalisationScheme","alphaMZ");
    break;
  case 3:
    SetParameter("EWRenormalisationScheme","alphaGF");
    break;
  default:
    THROW(fatal_error,"Unknown electroweak renormalisation scheme.");
    break;
  }

  // set couplings to unity (corrected in interface)
  GoSam_Interface::SetParameter("alpha", MODEL::aqed->Default());
//  GoSam_Interface::SetParameter("alpha_s", AlphaQCD());


  // set particle masses/widths
//  Flavour_Vector flavs(MODEL::s_model->IncludedFlavours());
//  for (size_t i(0); i<flavs.size(); ++i) {
//    SetParameter("mass("+ToString(flavs[i].Kfcode())+")", flavs[i].Mass());
//    SetParameter("width("+ToString(flavs[i].Kfcode())+")", flavs[i].Width());
//    if (flavs[i].IsFermion() && flavs[i].Mass()!=flavs[i].Yuk()) {
//      SetParameter("yukawa("+ToString(flavs[i].Kfcode())+")", flavs[i].Yuk());
//      if (flavs[i].IsQuark()) {
//        if (MODEL::s_model->ScalarNumber(std::string("YukawaScheme"))==1)
//          SetParameter("muy("+ToString(flavs[i].Kfcode())+")", Flavour(kf_h0).Mass(true));
//        else
//          SetParameter("muy("+ToString(flavs[i].Kfcode())+")", flavs[i].Yuk());
//      }
//    }
//  }

//  if (s_model->ComplexConstant("CKM_0_2")!=Complex(0.0,0.0) ||
//      s_model->ComplexConstant("CKM_2_0")!=Complex(0.0,0.0)) {
//    SetParameter("CKMorder", 3);
//  }
//  else if (s_model->ComplexConstant("CKM_1_2")!=Complex(0.0,0.0) ||
//	   s_model->ComplexConstant("CKM_2_1")!=Complex(0.0,0.0)) {
//    SetParameter("CKMorder", 2);
//  }
//  else if (s_model->ComplexConstant("CKM_0_1")!=Complex(0.0,0.0) ||
//	   s_model->ComplexConstant("CKM_1_0")!=Complex(0.0,0.0)) {
//    SetParameter("CKMorder", 1);
//  }
//  else {
//    SetParameter("CKMorder", 0);
//  }

  // set remaining GoSam parameters specified by user
  vector<string> parameters;
  reader.VectorFromFile(parameters,"GOSAM_PARAMETERS");
  for (size_t i=1; i<parameters.size(); i=i+2)
    SetParameter(parameters[i-1], parameters[i]);

  // instruct GoSam to return I instead of V
  if (s_vmode&2) SetParameter("VMode","I");

  char infostr[700];
  gosam_info(infostr);
  msg_Info()<<std::string(infostr)<<std::endl;

  MyStrStream cite;
  cite<<"The GoSam one-loop generator~\\cite{xxx:2010yy} "<<endl
      <<"has been used. "<<endl;
  rpa->gen.AddCitation(1,cite.str());

  return true;
}

int GoSam_Interface::RegisterProcess(const Subprocess_Info& is,
                                     const Subprocess_Info& fs)
{
  DEBUG_FUNC("");
  string shprocname(PHASIC::Process_Base::GenerateName(is,fs)),gsprocname("");
  Flavour_Vector isflavs(is.GetExternal());

  for (size_t i=0; i<isflavs.size(); ++i)
    gsprocname += ToString((long int)isflavs[i]) + " ";
  gsprocname += "-> ";
  Flavour_Vector fsflavs(fs.GetExternal());
  for (size_t i=0; i<fsflavs.size(); ++i)
    gsprocname += ToString((long int)fsflavs[i]) + " ";
  msg_Debugging()<<"looking for "<<shprocname<<" ("<<gsprocname<<")\n";

  int newlibs(0);
  int id(gosam_register_process(gsprocname.c_str(),&newlibs));

  if (newlibs) s_newlibs=true;

  if (s_procmap.find(id)==s_procmap.end())
    s_procmap[id]=shprocname;
  msg_Tracking()<<"GoSam_Interface process list:"<<std::endl;
  for (std::map<int,std::string>::const_iterator it=s_procmap.begin();
       it!=s_procmap.end();++it)
    msg_Tracking()<<it->first<<": "<<it->second<<std::endl;
  return id;
}

void GoSam_Interface::EvaluateLoop(int id, const Vec4D_Vector& moms,
                                   double& mu2,
                                   double& born, METOOLS::DivArrD& virt,
                                   double& accu)
{
  DEBUG_FUNC(id);
  double acc;
  Vec4D_Vector momenta(moms);
  Poincare cms(momenta[0]+momenta[1]);
  for (size_t i(0);i<momenta.size();++i) cms.Boost(momenta[i]);
  std::vector<double> pp(5*momenta.size());
  for (size_t i=0; i<momenta.size(); ++i) {
    pp[0+i*5]=momenta[i][0];
    pp[1+i*5]=momenta[i][1];
    pp[2+i*5]=momenta[i][2];
    pp[3+i*5]=momenta[i][3];
  }

  double mu(sqrt(mu2));
  std::vector<double> loop(4);
  gosam_evaluate_loop(&id, &pp[0], &mu, &loop[0], &accu);
  virt.Finite()=loop[2];
  virt.IR()=loop[1];
  virt.IR2()=loop[0];
  born=loop[3];
  msg_Debugging()<<"Born       = "<<born<<std::endl;
  msg_Debugging()<<"V_finite   = "<<virt.Finite()<<std::endl;
  msg_Debugging()<<"V_epsilon  = "<<virt.IR()<<std::endl;
  msg_Debugging()<<"V_epsilon2 = "<<virt.IR2()<<std::endl;
  msg_Debugging()<<"Accuracy   = "<<accu<<std::endl;
}

namespace GoSam {
  template <class ValueType>
  void HandleParameterStatus(int err, const std::string & key, ValueType value) {
    if (err==1) {
      msg_Debugging()<<"Setting GoSam parameter: "<<key<<" = "<<value<<endl;
    }
    else if (err==2) {
      msg_Info()<<"Unknown GoSam parameter: "<<key<<" = "<<ToString(value)
                <<std::endl;
    }
    else if (err==0) {
      std::string errorstring("Error setting GoSam parameter: "+key+" = "
                              +ToString(value));
      if (GoSam_Interface::ExitOnError()) THROW(fatal_error, errorstring)
      else                                msg_Error()<<errorstring<<std::endl;
    }
  }
}

void GoSam_Interface::SetParameter(const std::string & key, Complex value)
{
  int err(0);
  double re(value.real()), im(value.imag());
  gosam_setparameter_complex(key.c_str(), &re, &im, &err);
  GoSam::HandleParameterStatus(err,key,value);
}

void GoSam_Interface::SetParameter(const std::string & key, double value)
{
  int err(0);
  gosam_setparameter_double(key.c_str(), &value, &err);
  GoSam::HandleParameterStatus(err,key,value);
}

void GoSam_Interface::SetParameter(const std::string & key, int value)
{
  int err(0);
  gosam_setparameter_int(key.c_str(), &value, &err);
  GoSam::HandleParameterStatus(err,key,value);
}

void GoSam_Interface::SetParameter(const std::string & key, std::string value)
{
  int err(0);
  gosam_setparameter_string(key.c_str(), value.c_str(), &err);
  GoSam::HandleParameterStatus(err,key,value);
}

int GoSam_Interface::PerformTests()
{
  gosam_start();
  exh->AddTerminatorObject(this);
  return 1;
}

void GoSam_Interface::PrepareTerminate()
{
  gosam_finish();
}


DECLARE_GETTER(GoSam_Interface,"GoSam",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
                                  GoSam_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new GoSam::GoSam_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,GoSam_Interface>::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"Interface to the GoSam loop ME generator"; 
}

