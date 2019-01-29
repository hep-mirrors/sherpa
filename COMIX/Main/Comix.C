#ifndef COMIX_Main_Comix_H
#define COMIX_Main_Comix_H

#include "COMIX/Main/Process_Group.H"
#include "COMIX/Amplitude/Amplitude.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace MODEL { class Model_Base; }

namespace COMIX {

  class Single_Process;

  class Comix: public Process_Group, public PHASIC::ME_Generator_Base {

  private :

    std::vector<std::vector<Single_Process*> > m_umprocs;
    std::vector<PHASIC::Process_Base*>         m_rsprocs;

    time_t m_mets;

#ifdef USING__Threading
    CDBG_ME_TID_Vector m_cts;
#endif

    void RegisterDefaults() const;

    void PrintLogo(std::ostream &s);
    void PrintVertices();

  public :

    // constructor
    Comix();

    // destructor
    ~Comix();

    // member functions
    bool Initialize(MODEL::Model_Base *const model,
                    BEAM::Beam_Spectra_Handler *const beamhandler,
                    PDF::ISR_Handler *const isrhandler);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi,
                                            bool add);
    int PerformTests();
    bool NewLibraries();

  }; // end of class Comix

} // end of namespace COMIX

#endif

#include "COMIX/Main/Single_Process.H"
#include "COMIX/Main/Single_Dipole_Term.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "MODEL/Main/Model_Base.H"
//#include "REMNANTS/Main/Remnant_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/CXXFLAGS.H"

using namespace COMIX;
using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

Comix::Comix():
  ME_Generator_Base("Comix")
{
#ifdef USING__Threading
  p_cts=&m_cts;
#endif
}

Comix::~Comix() 
{
#ifdef USING__Threading
  for (size_t i(0);i<m_cts.size();++i) {
    CDBG_ME_TID *tid(m_cts[i]);
    tid->m_s=0;
    pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
    int tec(0);
    if ((tec=pthread_join(tid->m_id,NULL)))
      THROW(fatal_error,"Cannot join thread"+ToString(i));
    pthread_mutex_unlock(&tid->m_t_mtx);
    pthread_mutex_destroy(&tid->m_t_mtx);
    pthread_mutex_destroy(&tid->m_s_mtx);
    pthread_cond_destroy(&tid->m_t_cnd);
    pthread_cond_destroy(&tid->m_s_cnd);
  }
#endif
}

#define RED(ARG) om::red<<ARG<<om::reset
#define GREEN(ARG) om::green<<ARG<<om::reset
#define BLUE(ARG) om::blue<<ARG<<om::reset
#define YELLOW(ARG) om::brown<<ARG<<om::reset
#define BLACK(ARG) ARG

void Comix::PrintLogo(std::ostream &s)
{
  s<<"+----------------------------------+\n";
  s<<"|                                  |\n";
  s<<"|      "<<RED("CCC")<<"  "<<GREEN("OOO")<<"  "
   <<BLUE("M")<<"   "<<BLUE("M")<<" "<<BLACK("I")<<" "
   <<YELLOW("X")<<"   "<<YELLOW("X")<<"     |\n";
  s<<"|     "<<RED("C")<<"    "<<GREEN("O")<<"   "
   <<GREEN("O")<<" "<<BLUE("MM")<<" "<<BLUE("MM")
   <<" "<<BLACK("I")<<"  "<<YELLOW("X")<<" "
   <<YELLOW("X")<<"      |\n";
  s<<"|     "<<RED("C")<<"    "<<GREEN("O")
   <<"   "<<GREEN("O")<<" "<<BLUE("M")<<" "
   <<BLUE("M")<<" "<<BLUE("M")<<" "<<BLACK("I")
   <<"   "<<YELLOW("X")<<"       |\n";
  s<<"|     "<<RED("C")<<"    "<<GREEN("O")
   <<"   "<<GREEN("O")<<" "<<BLUE("M")<<"   "
   <<BLUE("M")<<" "<<BLACK("I")<<"  "
   <<YELLOW("X")<<" "<<YELLOW("X")<<"      |\n";
  s<<"|      "<<RED("CCC")<<"  "<<GREEN("OOO")
   <<"  "<<BLUE("M")<<"   "<<BLUE("M")<<" "
   <<BLACK("I")<<" "<<YELLOW("X")<<"   "
   <<YELLOW("X")<<"     |\n";
  s<<"|                                  |\n";
  s<<"+==================================+\n";
  s<<"|  Color dressed  Matrix Elements  |\n";
  s<<"|     http://comix.freacafe.de     |\n";
  s<<"|   please cite  JHEP12(2008)039   |\n";
  s<<"+----------------------------------+\n";
#ifdef USING__Threading
  s<<"Comix was compiled with thread support.\n";
#endif
  rpa->gen.AddCitation
    (1,"Comix is published under \\cite{Gleisberg:2008fv}.");
}

void Comix::PrintVertices()
{
  if (msg_LevelIsDebugging()) {
    msg_Out()<<METHOD<<"(): {\n\n   Implemented currents:\n\n";
    Current_Getter::PrintGetterInfo(msg_Out(),10);
    msg_Out()<<"\n   Implemented lorentz calculators:\n\n";
    LC_Getter::PrintGetterInfo(msg_Out(),10);
    msg_Out()<<"\n   Implemented color calculators:\n\n";
    CC_Getter::PrintGetterInfo(msg_Out(),10);
    msg_Out()<<"\n}\n";
  }
}

bool Comix::Initialize(MODEL::Model_Base *const model,
		       BEAM::Beam_Spectra_Handler *const beamhandler,
		       PDF::ISR_Handler *const isrhandler)
{
  RegisterDefaults();

  p_model=model;
  p_int->SetBeam(beamhandler); 
  p_int->SetISR(isrhandler);
  SetPSMasses();

  Scoped_Settings s{ Settings::GetMainSettings()["COMIX"] };
  s_partcommit = s["PARTIAL_COMMIT"].Get<int>();
  PrintLogo(msg->Info());
  PrintVertices();

  int helpi;

  helpi = s["VL_MODE"].Get<int>();
  Vertex::SetVLMode(helpi);

  double helpd;

#ifdef USING__Threading
  helpi = s["THREADS"].Get<int>();
  if (helpi>0) {
    m_cts.resize(helpi);
    for (size_t i(0);i<m_cts.size();++i) {
      CDBG_ME_TID *tid(new CDBG_ME_TID());
      m_cts[i] = tid;
      pthread_cond_init(&tid->m_s_cnd,NULL);
      pthread_cond_init(&tid->m_t_cnd,NULL);
      pthread_mutex_init(&tid->m_s_mtx,NULL);
      pthread_mutex_init(&tid->m_t_mtx,NULL);
      pthread_mutex_lock(&tid->m_t_mtx);
      tid->m_s=1;
      int tec(0);
      if ((tec=pthread_create(&tid->m_id,NULL,&Amplitude::TCalcJL,(void*)tid)))
	THROW(fatal_error,"Cannot create thread "+ToString(i));
    }
  }
#endif
  My_In_File::OpenDB
    (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
  return true;
}

void Comix::RegisterDefaults() const
{
  Scoped_Settings s{ Settings::GetMainSettings()["COMIX"] };
  s["PARTIAL_COMMIT"].SetDefault(0);
  s["PMODE"].SetDefault("D");
  s["WF_MODE"].SetDefault(0);  // wave function mode
  s["PG_MODE"].SetDefault(0);  // print graph mode
  s["VL_MODE"].SetDefault(0);  // vertex label mode
  s["N_GPL"].SetDefault(3);    // graphs per line
  s["THREADS"].SetDefault(0);  // number of threads
}

PHASIC::Process_Base *Comix::
InitializeProcess(const PHASIC::Process_Info &pi, bool add)
{
  if (p_model==NULL) return NULL;
  m_umprocs.push_back(std::vector<Single_Process*>());
  PHASIC::Process_Base *newxs(NULL);
  bool oneisgroup(pi.m_ii.IsGroup()||pi.m_fi.IsGroup());
  std::map<std::string,std::string> pmap;
  if (oneisgroup) {
    newxs = new Process_Group();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    newxs->Get<COMIX::Process_Base>()->SetCTS(p_cts);
    if (!newxs->Get<Process_Group>()->Initialize(&pmap,&m_umprocs.back())) {
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    newxs->Get<COMIX::Process_Base>()->SetGPath(pi.m_gpath);
    if (!newxs->Get<PHASIC::Process_Group>()->ConstructProcesses()) {
      My_In_File::CloseDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/",0);
      msg_Debugging()<<METHOD<<"(): Construct failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    My_In_File::CloseDB
      (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/",0);
    msg_Tracking()<<"Initialized '"<<newxs->Name()<<"'\n";
  }
  else {
    newxs = new Single_Process();
    newxs->SetGenerator(this);
    newxs->Init(pi,p_int->Beam(),p_int->ISR());
    newxs->Integrator()->SetHelicityScheme(pi.m_hls);
    newxs->Get<COMIX::Process_Base>()->SetModel(p_model);
    newxs->Get<COMIX::Process_Base>()->SetCTS(p_cts);
    newxs->Get<COMIX::Process_Base>()->SetGPath(pi.m_gpath);
    if (!newxs->Get<Single_Process>()->Initialize
	(&pmap,&m_umprocs.back(),m_blocks,m_nproc)) {
      My_In_File::CloseDB
	(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/",0);
      msg_Debugging()<<METHOD<<"(): Init failed for '"
		     <<newxs->Name()<<"'\n";
      delete newxs;
      return NULL;
    }
    if (!newxs->Get<Single_Process>()->MapProcess())
      if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
    My_In_File::CloseDB
      (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/",0);
  }
  if (add) Add(newxs,1);
  else m_rsprocs.push_back(newxs);
  return newxs;
}

int Comix::PerformTests()
{
  My_In_File::CloseDB
    (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
  if (!Tests()) return 0;
  for (size_t i=0;i<m_rsprocs.size();++i)
    if (!m_rsprocs[i]->Get<COMIX::Process_Base>()->Tests()) return false;
  return 1;
}

bool Comix::NewLibraries()
{
  return false;
}

DECLARE_GETTER(Comix,"Comix",ME_Generator_Base,ME_Generator_Key);

ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,Comix>::
operator()(const ME_Generator_Key &key) const
{
  return new Comix();
}

void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,Comix>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Comix ME generator"; 
}
