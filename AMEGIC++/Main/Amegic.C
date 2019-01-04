#ifndef AMEGIC_Main_Amegic_H
#define AMEGIC_Main_Amegic_H

#include "AMEGIC++/Main/Process_Group.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace AMEGIC {

  class Amegic: public Process_Group,
		public PHASIC::ME_Generator_Base {
  private :

    std::string  m_path, m_file;

    MODEL::Model_Base *p_mmodel;
    Amegic_Model      *p_amodel;

    std::vector<PHASIC::Process_Base*> m_rsprocs;

    void DrawLogo(std::ostream &ostr);

  public :

    // constructor
    Amegic();

    // destructor
    ~Amegic();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beamhandler,
		    PDF::ISR_Handler *const isrhandler);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi,
                                            bool add);
    int PerformTests();
    bool NewLibraries();

  };// end of class Amegic

}// end of namespace AMEGIC

#endif

#include "AMEGIC++/Main/Topology.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "MODEL/UFO/UFO_Model.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace ATOOLS;

void Amegic::DrawLogo(std::ostream &ostr)
{
  ostr<<"+-----------------------------------------+\n";
  ostr<<"|   X   X   X XXXX  XXX  XXX  XXX         |\n";
  ostr<<"|  X X  XX XX X    X      X  X     X   X  |\n";
  ostr<<"| X   X X X X XXX  X XXX  X  X    XXX XXX |\n";
  ostr<<"| XXXXX X   X X    X   X  X  X     X   X  |\n";
  ostr<<"| X   X X   X XXXX  XXX  XXX  XXX         |\n";
  ostr<<"+-----------------------------------------+\n";
  ostr<<"| please cite: JHEP 0202:044,2002         |\n";
  ostr<<"+-----------------------------------------+\n";
  rpa->gen.AddCitation
    (1,"Amegic is published under \\cite{Krauss:2001iv}.");
}

Amegic::Amegic():
  ME_Generator_Base("Amegic"), p_mmodel(NULL), p_amodel(NULL)
{
  DrawLogo(msg->Info());
  p_testmoms=NULL;
  p_gen=this;
}

Amegic::~Amegic()
{
  delete p_amodel;
}

bool Amegic::Initialize(const std::string &path,const std::string &file,
			MODEL::Model_Base *const model,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler)
{
  m_path=path;
  m_file=file;
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model) &&
      !reader.Get<int>("AMEGIC_ALLOW_UFO", 0))
    THROW(fatal_error, "AMEGIC can only be used in built-in models. Please use Comix for UFO models.");
  p_mmodel=model;
  p_amodel = new Amegic_Model(model);
  p_int->SetBeam(beamhandler);
  p_int->SetISR(isrhandler);
  SetPSMasses(&reader);


  double helpd; int helpi; std::string helps;
  size_t helpt; cs_itype::type helpcsit; sbt::subtype helpst;
  std::vector<int> helpvi;

  // N_color
  helpd = reader.Get("N_COLOR", 3.0,
                     "number of colours", METHOD);
  rpa->gen.SetVariable("N_COLOR", ToString(helpd));

  // dipole alpha, kappa, kt2 and gsplit
  helpd = reader.Get("DIPOLE_AMIN", Max(rpa->gen.Accu(),1.0e-8),
                     "dipole \\alpha_{cut}", METHOD);
  rpa->gen.SetVariable("DIPOLE_AMIN", ToString(helpd));

  helpd = reader.Get("DIPOLE_ALPHA", 1.0,
                     "dipole \\alpha_{max}", METHOD);
  rpa->gen.SetVariable("DIPOLE_ALPHA", ToString(helpd));

  helpd = reader.Get("DIPOLE_ALPHA_FF",
                     ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA")),
                     "FF dipole \\alpha_{max}", METHOD);
  rpa->gen.SetVariable("DIPOLE_ALPHA_FF", ToString(helpd));

  helpd = reader.Get("DIPOLE_ALPHA_FI",
                     ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA")),
                     "FI dipole \\alpha_{max}", METHOD);
  rpa->gen.SetVariable("DIPOLE_ALPHA_FI", ToString(helpd));

  helpd = reader.Get("DIPOLE_ALPHA_IF",
                     ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA")),
                     "IF dipole \\alpha_{max}", METHOD);
  rpa->gen.SetVariable("DIPOLE_ALPHA_IF", ToString(helpd));

  helpd = reader.Get("DIPOLE_ALPHA_II",
                     ToType<double>(rpa->gen.Variable("DIPOLE_ALPHA")),
                     "II dipole \\alpha_{max}", METHOD);
  rpa->gen.SetVariable("DIPOLE_ALPHA_II", ToString(helpd));

  helpd = reader.Get("DIPOLE_KAPPA", 2./3.,
                     "dipole \\kappa", METHOD);
  rpa->gen.SetVariable("DIPOLE_KAPPA", ToString(helpd));

  helpt = reader.Get("DIPOLE_NF_GSPLIT", Flavour(kf_jet).Size()/2,
                      "dipole N_f", METHOD);
  rpa->gen.SetVariable("DIPOLE_NF_GSPLIT", ToString(helpt));

  helpd = reader.Get("DIPOLE_KT2MAX", sqr(rpa->gen.Ecms()),
                     "dipole \\k_{T,max}^2", METHOD);
  rpa->gen.SetVariable("DIPOLE_KT2MAX", ToString(helpd));

  // parameters in dipole construction
  helpt = reader.Get("DIPOLE_COLLINEAR_VFF_SPLITTINGS", 1,
                     "collinear VFF splittings", METHOD);
  rpa->gen.SetVariable("DIPOLE_COLLINEAR_VFF_SPLITTINGS",ToString(helpt));

  helpt = reader.Get("DIPOLE_V_SUBTRACTION_MODE", 1,
                     "dipole V->VP subtraction mode", METHOD);
  rpa->gen.SetVariable("DIPOLE_V_SUBTRACTION_MODE",ToString(helpt));

  helpi = reader.Get("DIPOLE_PFF_IS_SPLIT_SCHEME", 1,
                     "split scheme in IS P->FF splittings", METHOD);
  rpa->gen.SetVariable("DIPOLE_PFF_IS_SPLIT_SCHEME",ToString(helpi));

  helpi = reader.Get("DIPOLE_PFF_FS_SPLIT_SCHEME", 0,
                     "split scheme in FS P->FF splittings", METHOD);
  rpa->gen.SetVariable("DIPOLE_PFF_FS_SPLIT_SCHEME",ToString(helpi));

  helpi = reader.Get("DIPOLE_PFF_IS_RECOIL_SCHEME", 0,
                     "recoil scheme in IS P->FF splittings", METHOD);
  rpa->gen.SetVariable("DIPOLE_PFF_IS_RECOIL_SCHEME",ToString(helpi));

  helpi = reader.Get("DIPOLE_PFF_FS_RECOIL_SCHEME", 0,
                     "recoil scheme in FS P->FF splittings", METHOD);
  rpa->gen.SetVariable("DIPOLE_PFF_FS_RECOIL_SCHEME",ToString(helpi));

  helpi = reader.Get("DIPOLE_IS_CLUSTER_TO_LEPTONS", 0,
                     "cluster dipoles to leptons ", METHOD);
  rpa->gen.SetVariable("DIPOLE_IS_CLUSTER_TO_LEPTONS",ToString(helpi));

  helpi = reader.Get("LIST_DIPOLES", 0,
                     "list dipoles", METHOD);
  rpa->gen.SetVariable("LIST_DIPOLES", ToString(helpi));

  // flavour restriction
  helps = reader.Get("DIPOLE_BORN_FLAVOUR_RESTRICTIONS", std::string(""),
                     "dipole underlying Born flavour restrictions", METHOD);
  rpa->gen.SetVariable("DIPOLE_BORN_FLAVOUR_RESTRICTIONS",helps);

  // nlo smearing parameters
  helpd = reader.Get("NLO_SMEAR_THRESHOLD", 0.,
                     "NLO smear threshold", METHOD);
  rpa->gen.SetVariable("NLO_SMEAR_THRESHOLD", ToString(helpd));

  helpd = reader.Get("NLO_SMEAR_POWER", 0.5,
                     "NLO smear power", METHOD);
  rpa->gen.SetVariable("NLO_SMEAR_POWER", ToString(helpd));

  // on-shell subtraction parameters
  helpi = reader.Get("DIPOLE_ONSHELL_SUBTRACTION", 0, "on-shell subtraction", METHOD);
  rpa->gen.SetVariable("DIPOLE_ONSHELL_SUBTRACTION",ToString(helpi));

  helpd = reader.Get("DIPOLE_ONSHELL_SUBTRACTION_WINDOW", 5.0, "on-shell subtraction window", METHOD);
  rpa->gen.SetVariable("DIPOLE_ONSHELL_SUBTRACTION_WINDOW",ToString(helpd));

  // OLP checks
  helpi = reader.Get("CHECK_BORN", 0,
                     "check Born against OLP", METHOD);
  rpa->gen.SetVariable("CHECK_BORN", ToString(helpi));

  helpi = reader.Get("CHECK_POLES", 0,
                     "check poles against OLP", METHOD);
  rpa->gen.SetVariable("CHECK_POLES", ToString(helpi));

  helpi = reader.Get("CHECK_FINITE", 0,
                     "check finite parts against OLP", METHOD);
  rpa->gen.SetVariable("CHECK_FINITE", ToString(helpi));

  helpi = reader.Get("CHECK_THRESHOLD", 0.0,
                     "threshold for checks", METHOD);
  rpa->gen.SetVariable("CHECK_THRESHOLD", ToString(helpi));

  // NLO options
  helpi =  reader.Get("LOOP_ME_INIT", 0,
                      "init the loop ME even when not needed", METHOD);
  rpa->gen.SetVariable("LOOP_ME_INIT", ToString(helpi));

  helpi = reader.Get("USR_WGT_MODE", 1,
                     "fill weight components for reweighting", METHOD);
  rpa->gen.SetVariable("USR_WGT_MODE", ToString(helpi));

  helpi = reader.Get("NLO_MUR_COEFFICIENT_FROM_VIRTUAL", 1,
                     "retrieve \\mu_R reweighting coefficients from OLP", METHOD);
  rpa->gen.SetVariable("NLO_MUR_COEFFICIENT_FROM_VIRTUAL", ToString(helpi));

  helpi = reader.Get("NLO_BVI_MODE", 0,
                     "BVI mode", METHOD);
  rpa->gen.SetVariable("NLO_BVI_MODE", ToString(helpi));

  helpcsit = reader.Get("NLO_IMODE", ToType<cs_itype::type>("IKP"),
                        "I-term component to calculate", METHOD);
  rpa->gen.SetVariable("NLO_IMODE", ToString(helpcsit));

  helpst = reader.Get("NLO_IPART", ToType<sbt::subtype>("QCD+QED"),
                      "I-term part", METHOD);
  rpa->gen.SetVariable("NLO_IPART", ToString(helpst));

  helpi = reader.Get("NLO_EPS_MODE", 0,
                     "prefactor mode", METHOD);
  rpa->gen.SetVariable("NLO_EPS_MODE", ToString(helpi));

  helpi = reader.Get("NLO_DR_MODE", 0,
                     "use dim. reduction if not dictated by OLP", METHOD);
  rpa->gen.SetVariable("NLO_DR_MODE", ToString(helpi));

  // general Amegic parameters
  helpi = reader.Get("AMEGIC_ALLOW_MAPPING", 1,
                     "allow process mapping", METHOD);
  rpa->gen.SetVariable("AMEGIC_ALLOW_MAPPING",ToString(helpi));

  helpi = reader.Get("AMEGIC_CHECK_LOOP_MAP", 0,
                     "check loop map", METHOD);
  rpa->gen.SetVariable("AMEGIC_CHECK_LOOP_MAP",ToString(helpi));

  helpi = reader.Get("AMEGIC_SORT_LOPROCESS", 1,
                     "sort LO processes", METHOD);
  rpa->gen.SetVariable("AMEGIC_SORT_LOPROCESS",ToString(helpi));

  helpi = reader.Get("AMEGIC_KEEP_ZERO_PROCS", 0,
                     "keep zero processes", METHOD);
  rpa->gen.Variable("AMEGIC_KEEP_ZERO_PROCS", ToString(helpi));

  helpi = reader.Get("AMEGIC_ME_LIBCHECK", 0,
                     "library check", METHOD);
  rpa->gen.SetVariable("AMEGIC_ME_LIBCHECK",ToString(helpi));

  helpi = reader.Get("AMEGIC_CUT_MASSIVE_VECTOR_PROPAGATORS", 1,
                     "cut massive vector propagators", METHOD);
  rpa->gen.SetVariable("AMEGIC_CUT_MASSIVE_VECTOR_PROPAGATORS",ToString(helpi));

  helpd = reader.Get("AMEGIC_DEFAULT_GAUGE", 1,
                     "gauge", METHOD);
  rpa->gen.SetVariable("AMEGIC_DEFAULT_GAUGE",ToString(helpd));

  AMEGIC::Process_Base::SetGauge(ToType<double>(rpa->gen.Variable("AMEGIC_DEFAULT_GAUGE")));

  s_partcommit = reader.Get("AMEGIC_PARTIAL_COMMIT", 0, "partial commit", METHOD);
  ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  return true;
}

PHASIC::Process_Base *Amegic::InitializeProcess(const PHASIC::Process_Info &pi,
                                                bool add)
{
  PHASIC::Process_Base *newxs(NULL);
  size_t nis(pi.m_ii.NExternal()), nfs(pi.m_fi.NExternal());
  std::string name(PHASIC::Process_Base::GenerateName(pi.m_ii,pi.m_fi));
  Topology top(nis+nfs);
  bool oneisgroup(pi.m_ii.IsGroup()||pi.m_fi.IsGroup());
  if (oneisgroup) {
    newxs = new AMEGIC::Process_Group();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    if (!newxs->Get<AMEGIC::Process_Group>()->
	InitAmplitude(p_amodel,&top)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (!newxs->Get<AMEGIC::Process_Group>()->ConstructProcesses()) {
      if (!s_partcommit)
	My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/",0);
      msg_Debugging()<<METHOD<<"(): Construct failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (!s_partcommit)
      My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/",0);
    newxs->Get<AMEGIC::Process_Group>()->WriteMappingFile();
    msg_Tracking()<<"Initialized '"<<newxs->Name()<<"'\n";
    if (msg_LevelIsTracking()) newxs->Get<AMEGIC::Process_Group>()->PrintProcessSummary();
  }
  else {
    newxs = GetProcess(pi);
    if (!newxs) return NULL;
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    p_testmoms = new Vec4D[newxs->NIn()+newxs->NOut()];
    if (!p_pinfo) {
      p_pinfo = Translate(pi);
      m_nin = newxs->NIn();
      m_flavs.clear();
      for (size_t i=0;i<m_nin;i++) 
	m_flavs.push_back(newxs->Flavours()[i]);
    }
    Phase_Space_Handler::TestPoint(p_testmoms,&newxs->Info(),this);
//    Vec4D sum;
//    Poincare lab(Vec4D(sqrt(10.0),0.0,0.0,1.0));
//    msg_Debugging()<<"After boost:\n";
//    for (size_t i(0);i<nis+nfs;++i) {
//      lab.Boost(p_testmoms[i]);
//      sum+=i<m_nin?-p_testmoms[i]:p_testmoms[i];
//      msg_Debugging()<<"  p["<<i<<"] = "<<p_testmoms[i]<<"\n";
//    }
//    msg_Debugging()<<"} -> sum = "<<sum<<"\n";
//    Poincare rot(Vec4D::ZVEC,Vec4D(sqrt(14.0),1.0,2.0,3.0));
//    msg_Debugging()<<"After rotation {\n";
//    for (size_t i(0);i<nis+nfs;++i) {
//      rot.Rotate(p_testmoms[i]);
//      sum+=i<m_nin?-p_testmoms[i]:p_testmoms[i];
//      msg_Debugging()<<"  p["<<i<<"] = "<<p_testmoms[i]<<"\n";
//    }
//    msg_Debugging()<<"} -> sum = "<<sum<<"\n";
    newxs->Get<AMEGIC::Process_Base>()->SetTestMoms(p_testmoms);
    newxs->Get<AMEGIC::Process_Base>()->SetPrintGraphs(pi.m_gpath);
    if (!newxs->Get<AMEGIC::Process_Base>()->
	InitAmplitude(p_amodel,&top,m_umprocs,m_errprocs)) {
      My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/",0);
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/",0);
  }
  if (add) Add(newxs,1);
  else m_rsprocs.push_back(newxs);
  newxs->SetGenerator(this);
  return newxs;
}

int Amegic::PerformTests()
{
  int tests(Process_Group::PerformTests());
  if (NewLibs()) return -1;
  for (size_t i(0);i<m_rsprocs.size();++i) 
    if (m_rsprocs[i]->Get<AMEGIC::Amegic_Base>()->NewLibs()) return -1;
  Minimize();
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");
  return tests;
}

bool Amegic::NewLibraries()
{
  if (NewLibs()) return true;
  for (size_t i(0);i<m_rsprocs.size();++i)
    if (m_rsprocs[i]->Get<AMEGIC::Amegic_Base>()->NewLibs()) return true;
  return false;
}

DECLARE_GETTER(Amegic,"Amegic",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,Amegic>::
operator()(const ME_Generator_Key &key) const
{
  return new Amegic();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,Amegic>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The AMEGIC++ ME generator"; 
}

