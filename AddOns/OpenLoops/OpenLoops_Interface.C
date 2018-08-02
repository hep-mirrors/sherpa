#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/UFO/UFO_Model.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include <algorithm>
#include <sys/stat.h>

#include "OpenLoops_Interface.H"

using namespace OpenLoops;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


extern "C" {

  void ol_welcome(char* str);
  void ol_set_init_error_fatal(int flag);
  int  ol_get_error();

  void ol_getparameter_double(const char* key, double val);
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
  void ol_evaluate_sc (int id, double* pp, int emitter, double* polvect, double* m2sc);
  void ol_evaluate_sc2(int id, double* pp, int emitter, double* polvect, double* m2sc);
  void ol_evaluate_cc (int id, double* pp, double* tree, double* m2cc, double *m2ewcc);
  void ol_evaluate_cc2(int id, double* pp, double* tree, double* m2cc, double *m2ewcc);
  void ol_evaluate_ccmatrix (int id, double* pp, double* tree, double* m2cc, double* m2ewcc);
  void ol_evaluate_ccmatrix2(int id, double* pp, double* tree, double* m2cc, double* m2ewcc);
}

// private static member definitions
std::string OpenLoops_Interface::s_olprefix = std::string("");
bool OpenLoops_Interface::s_ignore_model = false;
bool OpenLoops_Interface::s_exit_on_error = true;
std::vector<std::string> OpenLoops_Interface::s_evgen_params;

// private static member definitions
std::map<int,std::string> OpenLoops_Interface::s_procmap;
size_t OpenLoops_Interface::s_vmode;

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
  Default_Reader reader;
  s_ignore_model = reader.Get<int>("OL_IGNORE_MODEL", 0);
  s_exit_on_error = reader.Get<int>("OL_EXIT_ON_ERROR", 1);
  if (s_ignore_model) msg_Info()<<METHOD<<"(): OpenLoops will use the "
                                <<"Standard Model even if you set a "
                                <<"different model without warning."
                                <<std::endl;
  s_olprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/OpenLoops";
  if(stat(s_olprefix.c_str(),&st) != 0) s_olprefix = OPENLOOPS_PREFIX;
  s_olprefix = reader.Get<string>("OL_PREFIX", s_olprefix);
  msg_Info()<<"Initialising OpenLoops generator from "<<s_olprefix<<endl;
  s_vmode = reader.Get<int>("OL_VMODE", 0);
  msg_Tracking()<<METHOD<<"(): Set V-mode to "<<s_vmode<<endl;

  // load library dynamically
  s_loader->AddPath(s_olprefix+"/lib");
  s_loader->AddPath(s_olprefix+"/proclib");
  if (!s_loader->LoadLibrary("openloops"))
    THROW(fatal_error, "Failed to load libopenloops.");

  ol_set_init_error_fatal(0);

  // set OL verbosity
  std::string ol_verbosity = reader.Get<std::string>("OL_VERBOSITY", "0");
  SetParameter("verbose",ol_verbosity);

  // tell OL about the current model and check whether accepted
  if (!s_ignore_model) SetParameter("model", model->Name());

  // Propagate model parameters to OpenLoops
  if(dynamic_cast<UFO::UFO_Model*>(model))
    SetParametersUFO(model);
  else
    SetParametersSM(model);

  // set nf in alpha-s evolution
  int asnf0(isr->PDF(0)?isr->PDF(0)->ASInfo().m_nf:-1);
  if (asnf0==-1) asnf0=MODEL::as->Nf(1.e20);
  int asnf1(isr->PDF(1)?isr->PDF(1)->ASInfo().m_nf:-1);
  if (asnf1==-1) asnf1=MODEL::as->Nf(1.e20);
  if (asnf0==asnf1) SetParameter("minnf_alphasrun", asnf0);

  // set OL path
  SetParameter("install_path", s_olprefix.c_str());

#ifdef USING__MPI
    if (MPI::COMM_WORLD.Get_size()>1) SetParameter("splash","0");
#endif

  // set remaining OL parameters specified by user
  reader.ReadVector(s_evgen_params,"OL_PARAMETERS");
  for (size_t i=1; i<s_evgen_params.size(); i=i+2)
    SetParameter(s_evgen_params[i-1], s_evgen_params[i]);
  vector<string> int_params;
  reader.ReadVector(int_params,"OL_INTEGRATION_PARAMETERS");
  for (size_t i=1; i<int_params.size(); i=i+2)
    SetParameter(int_params[i-1], int_params[i]);

  if (s_vmode==2) SetParameter("ir_on",2);

  char welcomestr[700];
  ol_welcome(welcomestr);
  msg_Info()<<std::string(welcomestr)<<std::endl;

  MyStrStream cite;
  cite<<"The OpenLoops library~\\cite{Cascioli:2011va} of virtual"<<endl
      <<"matrix elements has been used. "<<endl;
  if (GetIntParameter("redlib1")==1 || GetIntParameter("redlib1")==7 ||
      GetIntParameter("redlib2")==1 || GetIntParameter("redlib2")==7) {
    cite<<"It is partly based on the tensor integral reduction described "<<endl
        <<"in~\\cite{Denner:2002ii,Denner:2005nn,Denner:2010tr}."<<endl;
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
        <<"in~\\cite{Mastrolia:2012bu,Peraro:2014cba}."<<endl;
  }
  rpa->gen.AddCitation(1,cite.str());
  return true;
}

// Propagate model parameters to OpenLoops in the standard model
void OpenLoops_Interface::SetParametersSM(const MODEL::Model_Base* model)
{
  // set ew scheme to as(mZ), irrespective of the Sherpa scheme,
  // we give parameters to OL as as(MZ) and masses
  SetParameter("ew_scheme",2);
  // ew-renorm-scheme to Gmu by default
  SetParameter("ew_renorm_scheme",1);

  // set particle masses/widths
  int tmparr[] = {kf_e, kf_mu, kf_tau, kf_u, kf_d, kf_s, kf_c, kf_b, kf_t,
                  kf_Wplus, kf_Z, kf_h0};
  vector<int> pdgids (tmparr, tmparr + sizeof(tmparr) / sizeof(tmparr[0]) );
  for (size_t i=0; i<pdgids.size(); ++i) {
    const int& id(pdgids[i]); const Flavour& flav(id);
    if (flav.Mass()>0.0) SetParameter("mass("+ToString(id)+")", flav.Mass());
    if (flav.Width()>0.0) SetParameter("width("+ToString(id)+")", flav.Width());
    if (flav.IsFermion() && flav.Yuk()>0.0 &&
        flav.Mass()!=flav.Yuk()) {
      SetParameter("yuk("+ToString(id)+")", flav.Yuk());
      if (flav.IsQuark()) { // not supported/needed for leptons
        if (model->ScalarNumber(std::string("YukawaScheme"))==1)
          SetParameter("muy("+ToString(id)+")", Flavour(kf_h0).Mass(true));
        else
          SetParameter("muy("+ToString(id)+")", flav.Yuk());
      }
    }
  }
  // Set CKM parameters
  if (model->ComplexConstant("CKM_0_2")!=Complex(0.0,0.0) ||
      model->ComplexConstant("CKM_2_0")!=Complex(0.0,0.0)) {
    SetParameter("ckmorder", 3);
  }
  else if (model->ComplexConstant("CKM_1_2")!=Complex(0.0,0.0) ||
      model->ComplexConstant("CKM_2_1")!=Complex(0.0,0.0)) {
    SetParameter("ckmorder", 2);
  }
  else if (model->ComplexConstant("CKM_0_1")!=Complex(0.0,0.0) ||
      model->ComplexConstant("CKM_1_0")!=Complex(0.0,0.0)) {
    SetParameter("ckmorder", 1);
  }
  else {
    SetParameter("ckmorder", 0);
  }
}

void OpenLoops_Interface::SetParametersUFO(const MODEL::Model_Base* model)
{
  // All external UFO parameters are stored in this map
  for(MODEL::ScalarConstantsMap::const_iterator it=model->ScalarConstants().begin();
      it!=model->ScalarConstants().end(); ++it)
    SetParameter(it->first, it->second);
}

void OpenLoops_Interface::SwitchMode(const int mode)
{
  for (size_t i=1; i<s_evgen_params.size(); i=i+2)
    SetParameter(s_evgen_params[i-1], s_evgen_params[i]);
}

int OpenLoops_Interface::RegisterProcess(const ATOOLS::Flavour_Vector& isflavs,
					 const ATOOLS::Flavour_Vector& fsflavs,
					 int amptype)
{
  PHASIC::Subprocess_Info ii;  PHASIC::Subprocess_Info fi;
  for (auto fl : isflavs)
    ii.m_ps.push_back(PHASIC::Subprocess_Info(fl));
  for (auto fl : fsflavs)
    fi.m_ps.push_back(PHASIC::Subprocess_Info(fl));
  return RegisterProcess(ii,fi,amptype);
}

int OpenLoops_Interface::RegisterProcess(const Subprocess_Info& is,
                                         const Subprocess_Info& fs,
                                         int amptype)
{
  DEBUG_FUNC("");
  string shprocname(PHASIC::Process_Base::GenerateName(is,fs)),olprocname("");
  Flavour_Vector isflavs(is.GetExternal());

  for (size_t i=0; i<isflavs.size(); ++i)
    olprocname += ToString((long int)isflavs[i]) + " ";
  olprocname += "-> ";
  Flavour_Vector fsflavs(fs.GetExternal());
  for (size_t i=0; i<fsflavs.size(); ++i)
    olprocname += ToString((long int)fsflavs[i]) + " ";
  msg_Debugging()<<"looking for "<<shprocname<<" ("<<olprocname<<")\n";

  int id(ol_register_process(olprocname.c_str(), amptype));
  if (s_procmap.find(id)==s_procmap.end())
    s_procmap[id]=shprocname;
  msg_Tracking()<<"OpenLoops_Interface process list:"<<std::endl;
  for (std::map<int,std::string>::const_iterator it=s_procmap.begin();
       it!=s_procmap.end();++it)
    msg_Tracking()<<it->first<<": "<<it->second<<std::endl;
  return id;
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

void OpenLoops_Interface::EvaluateLoop(int id, const Vec4D_Vector& momenta,
                                       double& res, METOOLS::DivArrD& virt)
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
  msg_Debugging()<<"Born       = "<<res<<std::endl;
  msg_Debugging()<<"V_finite   = "<<virt.Finite()<<std::endl;
  msg_Debugging()<<"V_epsilon  = "<<virt.IR()<<std::endl;
  msg_Debugging()<<"V_epsilon2 = "<<virt.IR2()<<std::endl;
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

double OpenLoops_Interface::EvaluateSpinCorrelator(int id, const Vec4D_Vector& momenta,
						   const Vec4D& polv,
						   size_t emitter, size_t spectator,
						   AmplitudeType type)
{
  vector<double> pp(5*momenta.size());
  for (size_t i=0; i<momenta.size(); ++i) {
    pp[0+i*5]=momenta[i][0];
    pp[1+i*5]=momenta[i][1];
    pp[2+i*5]=momenta[i][2];
    pp[3+i*5]=momenta[i][3];
  }
  double polvec[4] = {polv[0], polv[1], polv[2], polv[3]};
  std::vector<double> res(momenta.size(), 0.0);

  /* Add +1 to emitter index: fortran-style */
  if(type == Tree)
    ol_evaluate_sc(id, &pp[0], emitter+1,  polvec, &res[0]);
  else if(type == Loop2)
    ol_evaluate_sc2(id, &pp[0], emitter+1,  polvec, &res[0]);
  else
    THROW(fatal_error, "Unknown amplitude type");

  /* OpenLoops gets the sign wrong! */
  return -res[spectator];
}

double OpenLoops_Interface::EvaluateColorCorrelator(int id, const Vec4D_Vector& momenta,
						    size_t emitter, size_t spectator,
						    AmplitudeType type)
{
  vector<double> pp(5*momenta.size());
  for (size_t i=0; i<momenta.size(); ++i) {
    pp[0+i*5]=momenta[i][0];
    pp[1+i*5]=momenta[i][1];
    pp[2+i*5]=momenta[i][2];
    pp[3+i*5]=momenta[i][3];
  }
  size_t dim = momenta.size();
  vector<double> cc_qcd(dim*(dim-1)/2, 0.0);
  double cc_ew, dummy_tree;

  if(type == Tree)
    ol_evaluate_cc(id, &pp[0], &dummy_tree, &cc_qcd[0], &cc_ew);
  else if(type == Loop2)
    ol_evaluate_cc2(id, &pp[0], &dummy_tree, &cc_qcd[0], &cc_ew);
  else
    THROW(fatal_error, "Unknown amplitude type");

  if(emitter>spectator) std::swap(emitter, spectator);
  size_t index = emitter+spectator*(spectator-1)/2;

  return cc_qcd[index];
}

void OpenLoops_Interface::PopulateColorCorrelatorMatrix(int id, const Vec4D_Vector& momenta,
							double& born2, double* ccmatrix,
							AmplitudeType type)
{
  /* In principle no need to return the born2 (the squared
     uncorrelated ME) because it is proportional to the diagonal
     elements of the ccmatrix. Expose it nonetheless for validation
     purposes. */

  vector<double> pp(5*momenta.size());
  for (size_t i=0; i<momenta.size(); ++i) {
    pp[0+i*5]=momenta[i][0];
    pp[1+i*5]=momenta[i][1];
    pp[2+i*5]=momenta[i][2];
    pp[3+i*5]=momenta[i][3];
  }
  double dummy_cc_ew;

  if(type == Tree)
    ol_evaluate_ccmatrix(id, &pp[0], &born2, ccmatrix, &dummy_cc_ew);
  else if(type == Loop2)
    ol_evaluate_ccmatrix2(id, &pp[0], &born2, ccmatrix, &dummy_cc_ew);
  else
    THROW(fatal_error, "Unknown amplitude type");
}

double OpenLoops_Interface::GetDoubleParameter(const std::string & key)
{
  double value;
  ol_getparameter_double(key.c_str(), value);
  return value;
}
int OpenLoops_Interface::GetIntParameter(const std::string & key)
{
  int value;
  ol_getparameter_int(key.c_str(), &value);
  return value;
}
template <class ValueType>
void HandleParameterStatus(int err, const std::string & key, ValueType value)
{
  if (err==0) {
    msg_Debugging()<<"Setting OpenLoops parameter: "<<key<<" = "<<value<<endl;
  }
  else if (err==1) {
    std::string errorstring("Unknown OpenLoops parameter: "+key+" = "+ToString(value));
    if (OpenLoops_Interface::ExitOnError()) THROW(fatal_error, errorstring)
    else                                    msg_Error()<<errorstring<<std::endl;
  }
  else if (err==2) {
    std::string errorstring("Error setting OpenLoops parameter: "+key+" = "+ToString(value));
    if (OpenLoops_Interface::ExitOnError()) THROW(fatal_error, errorstring)
    else                                    msg_Error()<<errorstring<<std::endl;
  }
}
void OpenLoops_Interface::SetParameter(const std::string & key, double value)
{
  ol_setparameter_double(key.c_str(), &value);
  HandleParameterStatus(ol_get_error(), key, value);
}
void OpenLoops_Interface::SetParameter(const std::string & key, int value)
{
  ol_setparameter_int(key.c_str(), value);
  HandleParameterStatus(ol_get_error(), key, value);
}
void OpenLoops_Interface::SetParameter(const std::string & key, std::string value)
{
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

