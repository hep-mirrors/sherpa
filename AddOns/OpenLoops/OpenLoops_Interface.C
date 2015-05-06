#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>
#include <sys/stat.h>

#include "OpenLoops_Interface.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


extern "C" {

  void ol_welcome(char* str);
  void ol_set_init_error_fatal(int flag);
  int  ol_get_error();

  void ol_getparameter_double(const char* key, double* val);
  void ol_getparameter_int(const char* key, int* val);
  void ol_setparameter_double(const char* key, double* val);
  void ol_setparameter_int(const char* key, int val);
  void ol_setparameter_string(const char* key, const char* val);

  int ol_register_process(const char* process, int amptype);

  void ol_start();
  void ol_finish();

  void ol_evaluate_loop(int id, double* pp, double* m2l0, double* m2l1, double* acc);
  void ol_evaluate_tree(int id, double* pp, double* m2l0);
  void ol_evaluate_loop2(int id, double* pp, double* m2l0, double* acc);
}


namespace OpenLoops {

  std::string OpenLoops_Interface::s_olprefix;

  OpenLoops_Interface::~OpenLoops_Interface()
  {
    ol_finish();
  }

  bool OpenLoops_Interface::Initialize(const string &path,const string &file,
                                       MODEL::Model_Base *const model,
                                       BEAM::Beam_Spectra_Handler *const beam,
                                       PDF::ISR_Handler *const isr)
  {
    // find OL installation prefix with several overwrite options
    struct stat st;
    Data_Reader reader(" ",";","#","=");
    s_olprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/OpenLoops";
    if(stat(s_olprefix.c_str(),&st) != 0) s_olprefix = OPENLOOPS_PREFIX;
    s_olprefix = reader.GetValue<string>("OL_PREFIX", s_olprefix);
    msg_Info()<<"Initialising OpenLoops generator from "<<s_olprefix<<endl;

    // load library dynamically
    s_loader->AddPath(s_olprefix+"/lib");
    s_loader->AddPath(s_olprefix+"/proclib");
    if (!s_loader->LoadLibrary("openloops")) THROW(fatal_error, "Failed to load libopenloops.");

    ol_set_init_error_fatal(0);

    // tell OL about the current model
    SetParameter("model", MODEL::s_model->Name());

    // set particle masses/widths
    int tmparr[] = {kf_e, kf_mu, kf_tau, kf_u, kf_d, kf_s, kf_c, kf_b, kf_t, kf_Wplus, kf_Z, kf_h0};
    vector<int> pdgids (tmparr, tmparr + sizeof(tmparr) / sizeof(tmparr[0]) );
    for (size_t i=0; i<pdgids.size(); ++i) {
      if (Flavour(pdgids[i]).Mass()>0.0) SetParameter("mass("+ToString(pdgids[i])+")", Flavour(pdgids[i]).Mass());
      if (Flavour(pdgids[i]).Width()>0.0) SetParameter("width("+ToString(pdgids[i])+")", Flavour(pdgids[i]).Width());
    }

    if (s_model->ComplexConstant("CKM_0_2")!=Complex(0.0,0.0) ||
        s_model->ComplexConstant("CKM_2_0")!=Complex(0.0,0.0)) {
      SetParameter("ckmorder", 3);
    }
    else if (s_model->ComplexConstant("CKM_1_2")!=Complex(0.0,0.0) ||
        s_model->ComplexConstant("CKM_2_1")!=Complex(0.0,0.0)) {
      SetParameter("ckmorder", 2);
    }
    else if (s_model->ComplexConstant("CKM_0_1")!=Complex(0.0,0.0) ||
        s_model->ComplexConstant("CKM_1_0")!=Complex(0.0,0.0)) {
      SetParameter("ckmorder", 1);
    }
    else {
      SetParameter("ckmorder", 0);
    }
    SetParameter("install_path", s_olprefix.c_str());

    // set remaining OL parameters specified by user
    vector<string> parameters;
    reader.VectorFromFile(parameters,"OL_PARAMETERS");
    for (size_t i=1; i<parameters.size(); i=i+2) SetParameter(parameters[i-1], parameters[i]);

    char welcomestr[GetIntParameter("welcome_length")];
    ol_welcome(welcomestr);
    msg_Info()<<std::string(welcomestr)<<std::endl;

    MyStrStream cite;
    cite<<"The OpenLoops library~\\cite{Cascioli:2011va} of virtual"<<endl
        <<"matrix elements has been used. "<<endl;
    if (GetIntParameter("redlib1")==1 || GetIntParameter("redlib1")==7 ||
        GetIntParameter("redlib2")==1 || GetIntParameter("redlib2")==7) {
      cite<<"It is partly based on the tensor integral reduction "
          <<"described in~\\cite{Denner:2002ii,Denner:2005nn,"
          <<"Denner:2010tr,Denner:2014gla}."<<endl;
    }
    if (GetIntParameter("redlib1")==5 || GetIntParameter("redlib2")==5) {
      cite<<"It is partly based on the integrand reduction described "<<endl
          <<"in~\\cite{Ossola:2007ax,vanHameren:2010cp}."<<endl;
    }
    if (GetIntParameter("redlib1")==6 || GetIntParameter("redlib2")==6) {
      cite<<"It is partly based on the integrand reduction described "<<endl
          <<"in~\\cite{Mastrolia:2010nb,vanHameren:2010cp}."<<endl;
    }
    if (GetIntParameter("redlib1")==8 || GetIntParameter("redlib2")==8) {
      cite<<"It is partly based on the integrand reduction described "<<endl
         <<"in~\\cite{vanDeurzen:2013saa,Peraro:2014cba}."<<endl;
    }
    rpa->gen.AddCitation(1,cite.str());

    return true;
  }

  int OpenLoops_Interface::RegisterProcess(const Subprocess_Info& is,
                                           const Subprocess_Info& fs,
                                           int amptype)
  {
    string procname;

    Flavour_Vector isflavs(is.GetExternal());
    for (size_t i=0; i<isflavs.size(); ++i) procname += ToString((long int) isflavs[i]) + " ";

    procname += "-> ";

    Flavour_Vector fsflavs(fs.GetExternal());
    for (size_t i=0; i<fsflavs.size(); ++i) procname += ToString((long int) fsflavs[i]) + " ";

    return ol_register_process(procname.c_str(), amptype);
  }

  void OpenLoops_Interface::EvaluateTree(int id, const Vec4D_Vector& momenta, double& res)
  {
    vector<double> pp(5*momenta.size());
    for (size_t i=0; i<momenta.size(); ++i) {
      pp[0+i*5]=momenta[i][0];
      pp[1+i*5]=momenta[i][1];
      pp[2+i*5]=momenta[i][2];
      pp[3+i*5]=momenta[i][3];
    }

    ol_evaluate_tree(id, &pp[0], &res);
  }

  void OpenLoops_Interface::EvaluateLoop(int id, const Vec4D_Vector& momenta, double& res, METOOLS::DivArrD& virt)
  {
    double acc;
    vector<double> pp(5*momenta.size());
    for (size_t i=0; i<momenta.size(); ++i) {
      pp[0+i*5]=momenta[i][0];
      pp[1+i*5]=momenta[i][1];
      pp[2+i*5]=momenta[i][2];
      pp[3+i*5]=momenta[i][3];
    }

    vector<double> m2l1(3);
    ol_evaluate_loop(id, &pp[0], &res, &m2l1[0], &acc);
    virt.Finite()=m2l1[0];
    virt.IR()=m2l1[1];
    virt.IR2()=m2l1[2];
  }

  void OpenLoops_Interface::EvaluateLoop2(int id, const Vec4D_Vector& momenta, double& res)
  {
    double acc;
    vector<double> pp(5*momenta.size());
    for (size_t i=0; i<momenta.size(); ++i) {
      pp[0+i*5]=momenta[i][0];
      pp[1+i*5]=momenta[i][1];
      pp[2+i*5]=momenta[i][2];
      pp[3+i*5]=momenta[i][3];
    }

    ol_evaluate_loop2(id, &pp[0], &res, &acc);
  }


  double OpenLoops_Interface::GetDoubleParameter(const std::string & key) {
    double value;
    ol_getparameter_double(key.c_str(), &value);
    return value;
  }
  int OpenLoops_Interface::GetIntParameter(const std::string & key) {
    int value;
    ol_getparameter_int(key.c_str(), &value);
    return value;
  }
  template <class ValueType>
  void HandleParameterStatus(int err, const std::string & key, ValueType value) {
    if (err==0) {
      msg_Debugging()<<"Setting OpenLoops parameter: "<<key<<" = "<<value<<endl;
    }
    else if (err==1) {
      THROW(fatal_error, "Unknown OpenLoops parameter: "+key+" = "+ToString(value));
    }
    else if (err==2) {
      THROW(fatal_error, "Error setting OpenLoops parameter: "+key+" = "+ToString(value));
    }
  }
  void OpenLoops_Interface::SetParameter(const std::string & key, double value) {
    ol_setparameter_double(key.c_str(), &value);
    HandleParameterStatus(ol_get_error(), key, value);
  }
  void OpenLoops_Interface::SetParameter(const std::string & key, int value) {
    ol_setparameter_int(key.c_str(), value);
    HandleParameterStatus(ol_get_error(), key, value);
  }
  void OpenLoops_Interface::SetParameter(const std::string & key, std::string value) {
    ol_setparameter_string(key.c_str(), value.c_str());
    HandleParameterStatus(ol_get_error(), key, value);
  }


  int OpenLoops_Interface::PerformTests()
  {
    ol_start();
    exh->AddTerminatorObject(this);
    return 1;
  }

  void OpenLoops_Interface::PrepareTerminate()
  {
    ol_finish();
  }


}

using namespace OpenLoops;

DECLARE_GETTER(OpenLoops_Interface,"OpenLoops",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
                                  OpenLoops_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new OpenLoops::OpenLoops_Interface();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,OpenLoops_Interface>::
PrintInfo(ostream &str,const size_t width) const
{ 
  str<<"Interface to the OpenLoops loop ME generator"; 
}

