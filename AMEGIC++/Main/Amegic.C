#ifndef AMEGIC_Main_Amegic_H
#define AMEGIC_Main_Amegic_H

#include "AMEGIC++/Main/Process_Group.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace AMEGIC {

  class Amegic: public Process_Group, public PHASIC::ME_Generator_Base {

  private :

    MODEL::Model_Base *p_mmodel;
    Amegic_Model      *p_amodel;

    std::vector<PHASIC::Process_Base*> m_rsprocs;

    void RegisterDefaults() const;

    void DrawLogo(std::ostream &ostr);

  public :

    // constructor
    Amegic();

    // destructor
    ~Amegic();

    // member functions
    bool Initialize(MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beamhandler,
		    PDF::ISR_Handler *const isrhandler,
		    YFS::YFS_Handler *const yfshandler);
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
  ME_Generator_Base("Amegic"),
  p_mmodel(NULL), p_amodel(NULL)
{
  DrawLogo(msg->Info());
  p_testmoms=NULL;
  p_gen=this;
}

Amegic::~Amegic()
{
  delete p_amodel;
}

bool Amegic::Initialize(MODEL::Model_Base *const model,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler,
			YFS::YFS_Handler *const yfshandler)
{
  Settings& s = Settings::GetMainSettings();
  Scoped_Settings amegicsettings{ s["AMEGIC"] };
  RegisterDefaults();

  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)
      && !amegicsettings["ALLOW_UFO"].Get<int>()) {
    THROW(fatal_error, "AMEGIC can only be used in built-in models. Please use Comix for UFO models.");
  }

  p_mmodel=model;
  p_amodel = new Amegic_Model(model);
  p_int->SetBeam(beamhandler);
  p_int->SetISR(isrhandler);
  p_int->SetYFS(yfshandler);
  SetPSMasses();

  AMEGIC::Process_Base::SetGauge(amegicsettings["DEFAULT_GAUGE"].Get<int>());

  s_partcommit = amegicsettings["PARTIAL_COMMIT"].Get<int>();

  ATOOLS::MakeDir(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Amegic/");

  return true;
}

void Amegic::RegisterDefaults() const
{
  Scoped_Settings s{ Settings::GetMainSettings()["AMEGIC"] };
  s["ALLOW_UFO"].SetDefault(false);
  s["SORT_LOPROCESS"].SetDefault(true);
  s["ME_LIBCHECK"].SetDefault(false);
  s["CUT_MASSIVE_VECTOR_PROPAGATORS"].SetDefault(true);
  s["DEFAULT_GAUGE"].SetDefault(1);
  s["PARTIAL_COMMIT"].SetDefault(0);
  s["ALLOW_MAPPING"].SetDefault(1);
  s["CHECK_LOOP_MAP"].SetDefault(0);
  s["KEEP_ZERO_PROCS"].SetDefault(0);
  s["CHECK_BORN"].SetDefault(false);
  s["CHECK_POLES"].SetDefault(false);
  s["CHECK_FINITE"].SetDefault(false);
  s["CHECK_THRESHOLD"].SetDefault(0.0);
  s["LOOP_ME_INIT"].SetDefault(false);
  s["NLO_BVI_MODE"].SetDefault(0);
  s["NLO_EPS_MODE"].SetDefault(0);
  s["NLO_DR_MODE"].SetDefault(0);
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
    newxs->Init(pi,p_int->Beam(),p_int->ISR(),p_int->YFS());
    if (!newxs->Get<AMEGIC::Process_Group>()->
	InitAmplitude(p_amodel,&top)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (!newxs->Get<AMEGIC::Process_Group>()->ConstructProcesses()) {
      msg_Debugging()<<METHOD<<"(): Construct failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    newxs->Get<AMEGIC::Process_Group>()->WriteMappingFile();
    msg_Tracking()<<"Initialized '"<<newxs->Name()<<"'\n";
    if (msg_LevelIsTracking()) newxs->Get<AMEGIC::Process_Group>()->PrintProcessSummary();
  }
  else {
    newxs = GetProcess(pi);
    if (!newxs) return NULL;
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR(),p_int->YFS());
    p_testmoms = new Vec4D[newxs->NIn()+newxs->NOut()];
    if (!p_pinfo) {
      p_pinfo = Translate(pi);
      m_nin = newxs->NIn();
      m_flavs.clear();
      for (size_t i=0;i<m_nin;i++) 
	m_flavs.push_back(newxs->Flavours()[i]);
    }
    Phase_Space_Handler::TestPoint(p_testmoms,&newxs->Info(),this);
    PrepareTestMoms(p_testmoms,newxs->NIn(),newxs->NOut());
    newxs->Get<AMEGIC::Process_Base>()->SetTestMoms(p_testmoms);
    newxs->Get<AMEGIC::Process_Base>()->SetPrintGraphs(pi.m_gpath);
    if (!newxs->Get<AMEGIC::Process_Base>()->
	InitAmplitude(p_amodel,&top,m_umprocs,m_errprocs)) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
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

