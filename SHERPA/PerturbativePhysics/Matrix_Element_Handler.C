#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/RUsage.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Strings.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "PDF/Main/NLOMC_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

#include <unistd.h>
#include <cctype>

using namespace SHERPA;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

void Matrix_Element_Handler::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["OVERWEIGHT_THRESHOLD"].SetDefault(1e12);
  s["MEH_NLOADD"].SetDefault(1);
  s["MEH_EWADDMODE"].SetDefault(0);
  s["MEH_QCDADDMODE"].SetDefault(0);
  s["EVENT_SEED_MODE"].SetDefault(0);
  s["EVENT_SEED_FILE"].SetDefault(
                              "ran.stat." + rpa->gen.Variable("RNG_SEED"));
  s["EVENT_SEED_INCREMENT"].SetDefault(1);
  s["GENERATE_RESULT_DIRECTORY"].SetDefault(true);

  s["COLOUR_SCHEME"]
    .SetDefault(1)
    .SetReplacementList(cls::ColorSchemeTags());

  s["HELICITY_SCHEME"]
    .SetDefault(1)
    .SetReplacementList(hls::HelicitySchemeTags());

  s["NLO_IMODE"].SetDefault("IKP");

  s["PSI"]["ASYNC"].SetDefault(false);
}

void Matrix_Element_Handler::RegisterMainProcessDefaults(
    Scoped_Settings& procsettings)
{
  procsettings["Cut_Core"].SetDefault(0);
  procsettings["CKKW"].SetDefault("");
  procsettings.DeclareVectorSettingsWithEmptyDefault({
      "Decay", "DecayOS", "No_Decay" });
}

Matrix_Element_Handler::Matrix_Element_Handler(MODEL::Model_Base *model):
  m_gens(),
  p_proc(NULL), p_beam(NULL), p_isr(NULL), p_model(model),
  m_eventmode(0), m_hasnlo(0),
  p_shower(NULL), p_nlomc(NULL),
  m_sum(0.0),
  m_ranidx(0), m_fosettings(0), p_ranin(NULL), p_ranout(NULL),
  m_respath("Results"),
  m_seedmode(0),
  m_nloadd(1), m_ewaddmode(0), m_qcdaddmode(0),
  m_evtinfo(Weight_Info()),
  m_ovwth(10.), m_weightfactor(1.),
  m_nlomode(nlo_mode::none)
{
  Settings& s = Settings::GetMainSettings();
  RegisterDefaults();
  m_respath = s["RESULT_DIRECTORY"].Get<std::string>();
  m_respath=ShortenPathName(m_respath);
  if (m_respath[0]!='/' && s.GetPath()!="")
    m_respath=s.GetPath()+"/"+m_respath;
  m_eventmode = ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
  m_ovwth = s["OVERWEIGHT_THRESHOLD"].Get<double>();
  m_seedmode = s["EVENT_SEED_MODE"].Get<int>();
  m_nloadd = s["MEH_NLOADD"].Get<int>();
  m_ewaddmode = s["MEH_EWADDMODE"].Get<int>();
  m_qcdaddmode = s["MEH_QCDADDMODE"].Get<int>();
  std::string seedfile{ s["EVENT_SEED_FILE"].Get<std::string>() };
#ifdef USING__GZIP
  seedfile+=".gz";
#endif
  if (m_seedmode==1) {
#ifdef USING__GZIP
    p_ranin = new ATOOLS::igzstream(seedfile.c_str());
#else
    p_ranin = new std::ifstream(seedfile.c_str());
#endif
    if (p_ranin!=NULL && !p_ranin->good()) THROW
      (fatal_error,"Cannot initialize random generator status file");
  }
  else if (m_seedmode==2) {
#ifdef USING__GZIP
    p_ranout = new ATOOLS::ogzstream(seedfile.c_str());
#else
    p_ranout = new std::ofstream(seedfile.c_str());
#endif
    if (p_ranout!=NULL && !p_ranout->good()) THROW
      (fatal_error,"Cannot initialize random generator status file");
  }
  else if (m_seedmode==3) {
    const size_t incr{ s["EVENT_SEED_INCREMENT"].Get<size_t>() };
    ran->SetSeedStorageIncrement(incr);
  }
}

Matrix_Element_Handler::~Matrix_Element_Handler()
{
  if (p_ranin) delete p_ranin;
  if (p_ranout) delete p_ranout;
  for (size_t i=0;i<m_pmaps.size();++i) {
    for (std::map<ATOOLS::nlo_type::code,StringProcess_Map*>::const_iterator
	   pmit(m_pmaps[i]->begin());pmit!=m_pmaps[i]->end();++pmit)
      delete pmit->second;
    delete m_pmaps[i];
  }
  for (size_t i=0; i<m_procs.size(); ++i)
    if (dynamic_cast<MCatNLO_Process*>(m_procs[i])) delete m_procs[i];
  if (p_nlomc) delete p_nlomc;
}

void Matrix_Element_Handler::InitNLOMC()
{
  Settings& s = Settings::GetMainSettings();
  std::string nlomc((m_nlomode==nlo_mode::mcatnlo)?"MC@NLO":"");
  nlomc += "_" + Settings::GetMainSettings()["NLOMC_GENERATOR"].Get<std::string>();
  p_nlomc = NLOMC_Getter::GetObject(nlomc,NLOMC_Key(p_model,p_isr));
}

#ifdef USING__Threading
void *Matrix_Element_Handler::TCalculateTotalXSecs(void *arg)
{
  MEH_TID *tid((MEH_TID*)arg);
  pthread_cond_wait(&tid->m_s_cnd,&tid->m_s_mtx);
  tid->p_proc->SetLookUp(true);
  tid->p_proc->CalculateTotalXSec(tid->m_respath,false);
  tid->p_proc->SetLookUp(false);
  tid->p_proc->Integrator()->SetUpEnhance();
  tid->p_next->p_prev=tid->p_prev;
  tid->p_prev->p_next=tid->p_next;
  if (tid->p_next==tid) pthread_cond_signal(tid->p_t_cnd,tid->p_t_mtx);
  else pthread_cond_signal(&tid->p_next->m_s_cnd,&tid->p_next->m_s_mtx);
  return NULL;
}
#endif

bool Matrix_Element_Handler::CalculateTotalXSecs() 
{
  Settings& s = Settings::GetMainSettings();
#ifdef USING__Threading
  bool async = s["PSI"]["ASYNC"].Get<bool>();
#endif
  bool storeresults = s["GENERATE_RESULT_DIRECTORY"].Get<bool>();
  if (storeresults) {
    My_In_File::OpenDB(m_respath+"/");
  }
  bool okay(true);
#ifndef USING__Threading
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->SetLookUp(true);
    if (!m_procs[i]->CalculateTotalXSec(m_respath,false)) okay=false;
    m_procs[i]->SetLookUp(false);
    m_procs[i]->Integrator()->SetUpEnhance();
  }
#else
  pthread_mutex_t t_mtx;
  pthread_cond_t t_cnd;
  pthread_cond_init(&t_cnd,NULL);
  pthread_mutex_init(&t_mtx,NULL);
  pthread_mutex_lock(&t_mtx);
  MEH_TID_Vector cts;
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->AddMEHThread(cts,&Matrix_Element_Handler::TCalculateTotalXSecs);
  async = async && (cts.size() > 1);
  for (size_t i(0);i<cts.size();++i) {
    cts[i]->m_respath=m_respath;
    cts[i]->p_t_mtx=&t_mtx;
    cts[i]->p_t_cnd=&t_cnd;
    if (async) cts[i]->m_r=1;
  }
  pthread_cond_signal(&cts.front()->m_s_cnd,&cts.front()->m_s_mtx);
  pthread_cond_wait(&t_cnd,&t_mtx);
  for (size_t i(0);i<cts.size();++i) {
    MEH_TID *tid(cts[i]);
    if (pthread_join(tid->m_id,NULL))
      THROW(fatal_error,"Cannot join thread "+ToString(i));
    pthread_mutex_unlock(&tid->m_s_mtx);
    pthread_mutex_unlock(&tid->m_r_mtx);
    pthread_mutex_destroy(&tid->m_s_mtx);
    pthread_mutex_destroy(&tid->m_r_mtx);
    pthread_cond_destroy(&tid->m_s_cnd);
  }
  pthread_mutex_unlock(&t_mtx);
  pthread_mutex_destroy(&t_mtx);
  pthread_cond_destroy(&t_cnd);
#endif
  if (storeresults) My_In_File::CloseDB(m_respath+"/");
  return okay;
}

void Matrix_Element_Handler::SetRandomSeed()
{
  if (m_seedmode==1) {
    m_ranidx=ran->ReadInStatus(*p_ranin,m_ranidx);
    if (m_ranidx==std::string::npos) {
      msg_Error()<<METHOD<<"(): Status file read error. Abort."<<std::endl;
      Abort();
    }
  }
  else if (m_seedmode==2) {
    m_ranidx=ran->WriteOutStatus(*p_ranout,m_ranidx);
    if (m_ranidx==std::string::npos) {
      msg_Error()<<METHOD<<"(): Status file write error. Abort."<<std::endl;
      Abort();
    }
  }
}

bool Matrix_Element_Handler::GenerateOneEvent() 
{
  Return_Value::IncCall(METHOD);
  p_proc=NULL;
  if (m_seedmode!=3) SetRandomSeed();
  p_isr->SetPDFMember();

  // calculate total selection weight sum
  m_sum=0.0;
  for (size_t i(0);i<m_procs.size();++i)
    m_sum+=m_procs[i]->Integrator()->SelectionWeight(m_eventmode);

  // generate trial events until we accept one
  for (size_t n(1);true;++n) {
    rpa->gen.SetNumberOfTrials(rpa->gen.NumberOfTrials()+1);
    if (m_seedmode==3)
      ran->ResetToLastIncrementedSeed();
    if (!GenerateOneTrialEvent())
      continue;
    m_evtinfo.m_ntrial=n;
    return true;
  }
  return false;
}

bool Matrix_Element_Handler::GenerateOneTrialEvent()
{
  // select process
  double disc(m_sum*ran->Get()), csum(0.0);
  Process_Base *proc(NULL);
  for (size_t i(0);i<m_procs.size();++i) {
    if ((csum+=m_procs[i]->Integrator()->
          SelectionWeight(m_eventmode))>=disc) {
      proc=m_procs[i];
      break;
    }
  }

  // try to generate an event for the selected process
  if (proc==NULL) THROW(fatal_error,"No process selected");
  ATOOLS::Weight_Info *info=proc->OneEvent(m_eventmode);
  p_proc=proc->Selected();
  if (p_proc->Generator()==NULL)
    THROW(fatal_error,"No generator for process '"+p_proc->Name()+"'");
  if (p_proc->Generator()->MassMode()!=0)
    THROW(fatal_error,"Invalid mass mode. Check your PS interface.");
  if (info==NULL)
    return false;
  m_evtinfo=*info;
  delete info;

  // calculate weight factor and/or apply unweighting and weight threshold
  const auto sw = p_proc->Integrator()->SelectionWeight(m_eventmode) / m_sum;
  double enhance = p_proc->Integrator()->PSHandler()->Enhance();
  double wf(rpa->Picobarn()/sw/enhance);
  if (m_eventmode!=0) {
    const auto max = p_proc->Integrator()->Max();
    const auto disc = max * ran->Get();
    const auto abswgt = std::abs(m_evtinfo.m_weightsmap.Nominal());
    if (abswgt < disc) {
      return false;
    }
    if (abswgt > max * m_ovwth) {
      Return_Value::IncWarning(METHOD);
      msg_Info() << METHOD<<"(): Point for '" << p_proc->Name()
                 << "' exceeds maximum by "
                 << abswgt / max - 1.0 << "." << std::endl;
      m_weightfactor = m_ovwth;
      wf *= max * m_ovwth / abswgt;
    } else {
      m_weightfactor = abswgt / max;
      wf /= Min(1.0, m_weightfactor);
    }
  }

  // trial event is accepted, apply weight factor
  m_evtinfo.m_weightsmap["ME"]*=wf;
  if (p_proc->GetSubevtList()) {
    (*p_proc->GetSubevtList())*=wf;
    p_proc->GetSubevtList()->MultMEwgt(wf);
  }
  if (p_proc->GetMEwgtinfo()) (*p_proc->GetMEwgtinfo())*=wf;
  return true;
}

std::vector<Process_Base*> Matrix_Element_Handler::InitializeProcess(
    Process_Info pi, NLOTypeStringProcessMap_Map*& pmap)
{
  CheckInitialStateOrdering(pi);
  std::vector<Process_Base*> procs;
  std::set<Process_Info> initialized_pi_set;
  std::vector<Flavour_Vector> fls(pi.ExtractMPL());
  std::vector<int> fid(fls.size(),0);
  Flavour_Vector fl(fls.size());
  for (size_t i(0);i<fid.size();++i) fl[i]=fls[i][0];
  for (size_t hc(fid.size()-1);fid[0]<fls[0].size();) {
    if(fid[hc]==fls[hc].size()){fid[hc--]=0;++fid[hc];continue;}
    fl[hc]=fls[hc][fid[hc]];if(hc<fid.size()-1){++hc;continue;}
    Flavour_Vector cfl(fl);
    size_t n{0};
    pi.m_ii.SetExternal(cfl, n);
    pi.m_fi.SetExternal(cfl, n);
    Process_Base::SortFlavours(pi,1);
    if (initialized_pi_set.find(pi)==initialized_pi_set.end()) {
      initialized_pi_set.insert(pi);
      std::vector<Process_Base*> cp=InitializeSingleProcess(pi,pmap);
      procs.insert(procs.end(),cp.begin(),cp.end());
    }
    ++fid[hc];
  }
  return procs;
}

std::vector<Process_Base*> Matrix_Element_Handler::InitializeSingleProcess
(const Process_Info &pi,NLOTypeStringProcessMap_Map *&pmap)
{
  std::vector<Process_Base*> procs;
  if (pi.m_fi.NLOType()==nlo_type::lo) {
    Process_Base *proc(m_gens.InitializeProcess(pi, true));
    if (proc) {
      m_procs.push_back(proc);
      procs.push_back(proc);
      if (pmap==NULL) {
	m_pmaps.push_back(new NLOTypeStringProcessMap_Map());
	pmap=m_pmaps.back();
      }
      m_procs.back()->FillProcessMap(pmap);
    }
    return procs;
  }
  else {
    if (m_nlomode==nlo_mode::mcatnlo) {
      m_hasnlo=3;
      if (p_nlomc==NULL) InitNLOMC();
      if (pmap==NULL) {
	m_pmaps.push_back(new NLOTypeStringProcessMap_Map());
	pmap=m_pmaps.back();
      }
      MCatNLO_Process *proc=new MCatNLO_Process(m_gens,pmap);
      proc->Init(pi,p_beam,p_isr);
      if ((*proc)[0]==NULL) {
	delete proc;
	return procs;
      }
      if (!p_shower->GetShower())
        THROW(fatal_error,"Shower needs to be set for MC@NLO");
      proc->SetShower(p_shower->GetShower());
      proc->SetNLOMC(p_nlomc);
      m_procs.push_back(proc);
      procs.push_back(proc);
      return procs;
    }
    else if (m_nlomode==nlo_mode::fixedorder) {
      m_hasnlo=1;
      if (pi.m_fi.NLOType()&(nlo_type::vsub|nlo_type::loop|nlo_type::born)) {
	Process_Info rpi(pi);
	rpi.m_fi.SetNLOType(pi.m_fi.NLOType()&(nlo_type::vsub|nlo_type::loop|
					       nlo_type::born));
	if (m_nloadd) {
	  if (rpi.m_fi.m_nlocpl.size()<2) THROW(fatal_error,"NLO_Order not set.");
	  for (int i(0);i<2;++i) {
	    rpi.m_maxcpl[i]+=rpi.m_fi.m_nlocpl[i];
	    rpi.m_mincpl[i]+=rpi.m_fi.m_nlocpl[i];
	  }
	}
	procs.push_back(m_gens.InitializeProcess(rpi,true));
	if (procs.back()==NULL) {
	  msg_Error()<<"No such process:\n"<<rpi<<std::endl;
	  THROW(critical_error,"Failed to intialize process");
	}
      }
      if (pi.m_fi.NLOType()&nlo_type::real ||
          pi.m_fi.NLOType()&nlo_type::rsub) {
        // if real or rsub is requested, the extra jet is not yet contained
        // in the process info, but has to be added here
        Process_Info rpi(pi);
	rpi.m_fi.SetNLOType(pi.m_fi.NLOType()&(nlo_type::real|nlo_type::rsub));
	rpi.m_integrator=rpi.m_rsintegrator;
	rpi.m_megenerator=rpi.m_rsmegenerator;
	rpi.m_itmin=rpi.m_rsitmin;
	if (m_nloadd) {
	  if (rpi.m_fi.m_nlocpl.size()<2) THROW(fatal_error,"NLO_Order not set.");
	  for (int i(0);i<2;++i) {
	    rpi.m_maxcpl[i]+=rpi.m_fi.m_nlocpl[i];
	    rpi.m_mincpl[i]+=rpi.m_fi.m_nlocpl[i];
	  }
	  if (pi.m_fi.m_nlocpl[0]==0. && pi.m_fi.m_nlocpl[1]==1.) {
	    if (m_ewaddmode==0)
	      rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_ewjet,"",""));
	    else if (m_ewaddmode==1)
	      rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
	    else if (m_ewaddmode==2)
	      rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
	    else
	      THROW(fatal_error,"Unknown MEH_EWADDMODE.");
	  }
	  else if (pi.m_fi.m_nlocpl[0]==1. && pi.m_fi.m_nlocpl[1]==0.) {
	    if (m_qcdaddmode==0)
	      rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_ewjet,"",""));
	    else if (m_qcdaddmode==1)
	      rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_gluon,"",""));
	    else if (m_qcdaddmode==2)
	      rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
	    else
	      THROW(fatal_error,"Unknown MEH_QCDADDMODE.");
	  }
	  else THROW(not_implemented,"Cannot do NLO QCD+EW yet.");
	}
        procs.push_back(m_gens.InitializeProcess(rpi,true));
	if (procs.back()==NULL) {
	  msg_Error()<<"No such process:\n"<<rpi<<std::endl;
	  THROW(critical_error,"Failed to intialize process");
	}
      }
      if (pmap==NULL) {
	m_pmaps.push_back(new NLOTypeStringProcessMap_Map());
	pmap=m_pmaps.back();
      }
      for (size_t i(0);i<procs.size();i++) {
	m_procs.push_back(procs[i]);
	m_procs.back()->FillProcessMap(pmap);
      }
      if (m_fosettings==0) {
        Settings& s = Settings::GetMainSettings();
	m_fosettings=1;
	if (p_shower->GetShower())
	  p_shower->GetShower()->SetOn(false);
        s["FRAGMENTATION"].OverrideScalar<std::string>("None");
        s["MI_HANDLER"].OverrideScalar<std::string>("None");
        if (!s["MI_HANDLER"].IsCustomised()) {
          s["BEAM_REMNANTS"].OverrideScalar<bool>(false);
        } else {
          Scoped_Settings kperp_settings{ s["INTRINSIC_KPERP"] };
          kperp_settings["MEAN"].OverrideScalar<double>(0.0);
          kperp_settings["SIGMA"].OverrideScalar<double>(0.0);
	}
        Scoped_Settings meqedsettings{ s["ME_QED"] };
        if (!meqedsettings["ENABLED"].IsCustomised()) {
          meqedsettings["ENABLED"].OverrideScalar<bool>(false);
        }
      }
    }
    else THROW(fatal_error,"NLO_Mode "+ToString(m_nlomode)+" unknown.");
  }
  return procs;
}

void Matrix_Element_Handler::CheckInitialStateOrdering(const Process_Info& pi)
{
  auto cpi = pi;
  Process_Base::SortFlavours(cpi, 1);
  if (cpi.m_ii == pi.m_ii) {
  } else {
    msg_Error() << ATOOLS::om::red << "\n\nERROR:" << ATOOLS::om::reset
      << " Wrong ordering of initial-state particles detected.\n"
      << "Please re-order the initial state in your Process definition(s) "
      << "like this:\n  ";
    cpi.m_ii.PrintFlavours(msg_Error());
    msg_Error() << " ->  ";
    pi.m_fi.PrintFlavours(msg_Error());
    msg_Error() << "\nYou may need to adjust your other beam-specific "
      << "parameters accordingly.\n";
    exit(-1);
  }
}

int Matrix_Element_Handler::InitializeProcesses(
  BEAM::Beam_Spectra_Handler* beam, PDF::ISR_Handler* isr)
{
  p_beam=beam;
  p_isr=isr;
  if (!m_gens.InitializeGenerators(p_model,beam,isr)) return false;
  double rbtime(ATOOLS::rpa->gen.Timer().RealTime());
  double btime(ATOOLS::rpa->gen.Timer().UserTime());
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  BuildProcesses();
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  if (msg_LevelIsTracking()) msg_Info()<<"Process initialization";
  double retime(ATOOLS::rpa->gen.Timer().RealTime());
  double etime(ATOOLS::rpa->gen.Timer().UserTime());
  size_t rss(GetCurrentRSS());
  msg_Info()<<" done ( "<<rss/(1<<20)<<" MB, "
	    <<FormatTime(size_t(retime-rbtime))<<" / "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
  if (m_procs.empty() && m_gens.size()>0)
    THROW(normal_exit,"No hard process found");
  msg_Info()<<METHOD<<"(): Performing tests "<<std::flush;
  rbtime=retime;
  btime=etime;
  int res(m_gens.PerformTests());
  retime=ATOOLS::rpa->gen.Timer().RealTime();
  etime=ATOOLS::rpa->gen.Timer().UserTime();
  rss=GetCurrentRSS();
  msg_Info()<<" done ( "<<rss/(1<<20)<<" MB, "
	    <<FormatTime(size_t(retime-rbtime))<<" / "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
  msg_Debugging()<<METHOD<<"(): Processes {\n";
  msg_Debugging()<<"  m_procs:\n";
  for (size_t i(0);i<m_procs.size();++i) 
    msg_Debugging()<<"    "<<m_procs[i]->Name()<<" -> "<<m_procs[i]<<"\n";
  msg_Debugging()<<"}\n";
  msg_Info()<<METHOD<<"(): Initializing scales "<<std::flush;
  rbtime=retime;
  btime=etime;
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  for (size_t i=0; i<m_procs.size(); ++i) m_procs[i]->InitScale();
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  retime=ATOOLS::rpa->gen.Timer().RealTime();
  etime=ATOOLS::rpa->gen.Timer().UserTime();
  rss=GetCurrentRSS();
  msg_Info()<<" done ( "<<rss/(1<<20)<<" MB, "
	    <<FormatTime(size_t(retime-rbtime))<<" / "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
  if (m_gens.NewLibraries())
    THROW(normal_exit,"Source code created. Run './makelibs' to compile.");
  return res;
}

void Matrix_Element_Handler::BuildProcesses()
{
  Settings& s = Settings::GetMainSettings();

  // init processes
  msg_Info()<<METHOD<<"(): Looking for processes "<<std::flush;
  if (msg_LevelIsTracking()) msg_Info()<<"\n";
  if (!m_gens.empty() && s["PROCESSES"].GetItemsCount() == 0) {
    if (!msg_LevelIsTracking()) msg_Info()<<"\n";
      THROW(missing_input, std::string{"Missing PROCESSES definition.\n\n"} +
                               Strings::ProcessesSyntaxExamples);
  }

  // iterate over processes in the settings
  for (auto& proc : s["PROCESSES"].GetItems()) {
    const auto keys = proc.GetKeys();
    if (keys.size() != 1) {
      if (!msg_LevelIsTracking()) msg_Info()<<"\n";
      THROW(invalid_input, std::string{"Invalid PROCESSES definition.\n\n"} +
                               Strings::ProcessesSyntaxExamples);
    }
    const std::string& name = keys[0];
    auto procsettings = proc[name];
    RegisterMainProcessDefaults(procsettings);
    Single_Process_List_Args args;
    ReadFinalStateMultiIndependentProcessSettings(name, procsettings, args);
    ReadFinalStateMultiSpecificProcessSettings(procsettings, args);
    BuildSingleProcessList(args);
    if (msg_LevelIsDebugging()) {
      msg_Indentation(4);
      msg_Out()<<m_procs.size()<<" process(es) found ..."<<std::endl;
      for (unsigned int i=0; i<m_procs.size(); ++i) {
        msg_Out()<<m_procs[i]->Name();
        if (m_procs[i]->IsGroup())
          msg_Out()<<" has subprocesses ...";
        msg_Out()<<std::endl;
      }
    }
  }
}

void Matrix_Element_Handler::ReadFinalStateMultiIndependentProcessSettings(
  const std::string& procname, Scoped_Settings proc, Single_Process_List_Args& args)
{
  // fill process info
  Settings& s = Settings::GetMainSettings();
  args.pi.m_scale = s["SCALES"].Get<std::string>();
  const auto couplings = s["COUPLINGS"].GetVector<std::string>();
  args.pi.m_coupling = MakeString(couplings, 0);
  args.pi.m_kfactor = s["KFACTOR"].Get<std::string>();
  args.pi.m_cls = (cls::scheme)s["COLOUR_SCHEME"].Get<int>();
  args.pi.m_hls = (hls::scheme)s["HELICITY_SCHEME"].Get<int>();
  std::vector<long int> nodecaykfcs{ proc["No_Decay"].GetVector<long int>() };
  for (const auto& kfc : nodecaykfcs)
    args.pi.m_nodecs.push_back(Flavour(std::abs(kfc),kfc<0));
  args.pi.p_gens=&m_gens;

  // fill initial/final state
  size_t pos(procname.find("->"));
  if (pos==std::string::npos)
    THROW(fatal_error, "Process name must be of the form `a b -> x y ...'.");
  args.ini = procname.substr(0,pos);
  args.fin = procname.substr(pos+2);

  // fill decay tags
  std::vector<std::string> helpsv;
  helpsv = proc["Decay"].GetVector<std::string>();
  args.dectags.insert(args.dectags.end(), helpsv.begin(), helpsv.end());
  helpsv = proc["DecayOS"].GetVector<std::string>();
  for (const auto& tag : helpsv)
    args.dectags.push_back("Z" + tag);

  // build process block info
  args.pbi.m_cutcore = proc["Cut_Core"].Get<int>();
  args.pbi.m_selectors = proc["Selectors"];

  // make modifications for multi-jet merging
  std::string ckkw{ proc["CKKW"].Get<std::string>() };
  if (ckkw != "") {
    if (p_shower==NULL || p_shower->GetShower()==NULL)
      THROW(fatal_error,"Invalid shower generator");
    args.pi.m_ckkw=1;
    args.pbi.m_gycut=ckkw;
  }
}

void Matrix_Element_Handler::ReadFinalStateMultiSpecificProcessSettings(
  Scoped_Settings proc,
  Single_Process_List_Args& args,
  std::string rawrange) const
{
  // iterate over settings in each process that can be either given for all
  // final-state multiplicities or only for some final-state multiplicities;
  // top-level items are applied to all final-state multiplicities (indicated
  // by the wildcard passed to ReadProcessSettings, whereas items that are
  // scoped under an integer key (e.g. "3: { ... }") or a continuous integer
  // range (e.g.  "4-7: { ... }") only apply to the corresponding final-state
  // multiplicities; ReadProcessSettings will recursively read all top-level
  // and scoped settings

  // parse final-state multiplicity range, allowed syntaxes are:
  // - single final-state multiplicity, e.g. "2"
  // - final-state multiplicity range, e.g. "2-4"
  // - one of the above, but prefixed with "x->", where x is any integer;
  //   this can be used to be more explicit, i.e. by writing "2->2-4", but note
  //   that the number before "->" is completely ignored by the parsing below
  // finally, the wildcard "-" stands for all multis, for which we leave
  // begin=0 and end=\infty
  auto range = std::make_pair(0, std::numeric_limits<size_t>::max());
  if (rawrange != "-") {
    const auto ranglepos = rawrange.find('>');
    if (ranglepos != std::string::npos)
      rawrange = rawrange.substr(ranglepos + 1);
    const auto delimiterpos = rawrange.find('-');
    range.first = std::stoul(rawrange.substr(0, delimiterpos));
    if (delimiterpos == std::string::npos)
      range.second = range.first;
    else
      range.second = std::stoul(rawrange.substr(delimiterpos + 1));
  }

  const auto nf = proc.GetIndex();
  for (auto rawsubkey : proc.GetKeys()) {

    // resolve tags which we allow in final-state multiplicity keys
    auto subkey = rawsubkey;
    proc.ReplaceTags(subkey);

    // recurse if a multiplicity specification is encountered
    if (std::isdigit(subkey[0])) {
      ReadFinalStateMultiSpecificProcessSettings(proc[rawsubkey], args, subkey);
      continue;
    }

    // ignore certain settings that are not to be handled by ExtractMPvalues
    // below
    if (subkey == "Selectors"
        || subkey == "Cut_Core"
        || subkey == "Decay"
        || subkey == "DecayOS"
        || subkey == "No_Decay"
        || subkey == "CKKW")
      continue;

    // read value (and potentially do some pre-processing for non-scalar
    // settings)
    std::string value;
    if (subkey == "Order"
        || subkey == "Max_Order"
        || subkey == "Min_Order"
        || subkey == "Amplitude_Order"
        || subkey == "Max_Amplitude_Order"
        || subkey == "Min_Amplitude_Order"
        || subkey == "NLO_Order") {
      // translate back into a single string to use ExtractMPvalues below
      value = MakeOrderString(proc[rawsubkey]);
    } else {
      value = proc[rawsubkey].SetDefault("").Get<std::string>();
    }

    if (subkey == "Order")                   ExtractMPvalues(value, range, nf, args.pbi.m_vcpl);
    else if (subkey == "Max_Order")          ExtractMPvalues(value, range, nf, args.pbi.m_vmaxcpl);
    else if (subkey == "Min_Order")          ExtractMPvalues(value, range, nf, args.pbi.m_vmincpl);
    else if (subkey == "Amplitude_Order")    ExtractMPvalues(value, range, nf, args.pbi.m_vacpl);
    else if (subkey == "Max_Amplitude_Order") ExtractMPvalues(value, range, nf, args.pbi.m_vmaxacpl);
    else if (subkey == "Min_Amplitude_Order") ExtractMPvalues(value, range, nf, args.pbi.m_vminacpl);
    else if (subkey == "Scales")             ExtractMPvalues(value, range, nf, args.pbi.m_vscale);
    else if (subkey == "Couplings")          ExtractMPvalues(value, range, nf, args.pbi.m_vcoupl);
    else if (subkey == "KFactor")            ExtractMPvalues(value, range, nf, args.pbi.m_vkfac);
    else if (subkey == "Y_Cut")              ExtractMPvalues(value, range, nf, args.pbi.m_vycut);
    else if (subkey == "Min_N_Quarks")       ExtractMPvalues(value, range, nf, args.pbi.m_vnminq);
    else if (subkey == "Max_N_Quarks")       ExtractMPvalues(value, range, nf, args.pbi.m_vnmaxq);
    else if (subkey == "Print_Graphs")       ExtractMPvalues(value, range, nf, args.pbi.m_vgpath);
    else if (subkey == "Name_Suffix")        ExtractMPvalues(value, range, nf, args.pbi.m_vaddname);
    else if (subkey == "Special")            ExtractMPvalues(value, range, nf, args.pbi.m_vspecial);
    else if (subkey == "Enable_MHV")         ExtractMPvalues(value, range, nf, args.pbi.m_vamegicmhv);
    else if (subkey == "Min_N_TChannels")    ExtractMPvalues(value, range, nf, args.pbi.m_vntchan);
    else if (subkey == "Max_N_TChannels")    ExtractMPvalues(value, range, nf, args.pbi.m_vmtchan);
    else if (subkey == "Integration_Error")  ExtractMPvalues(value, range, nf, args.pbi.m_vmaxerr);
    else if (subkey == "Max_Epsilon")        ExtractMPvalues(value, range, nf, args.pbi.m_vmaxeps);
    else if (subkey == "RS_Enhance_Factor")  ExtractMPvalues(value, range, nf, args.pbi.m_vrsefac);
    else if (subkey == "Enhance_Factor")     ExtractMPvalues(value, range, nf, args.pbi.m_vefac);
    else if (subkey == "Enhance_Function")   ExtractMPvalues(value, range, nf, args.pbi.m_vefunc);
    else if (subkey == "Enhance_Observable") ExtractMPvalues(value, range, nf, args.pbi.m_veobs);
    else if (subkey == "NLO_Mode")           ExtractMPvalues(value, range, nf, args.pbi.m_vnlomode);
    else if (subkey == "NLO_Part")           ExtractMPvalues(value, range, nf, args.pbi.m_vnlopart);
    else if (subkey == "NLO_Order")          ExtractMPvalues(value, range, nf, args.pbi.m_vnlocpl);
    else if (subkey == "Subdivide_Virtual")  ExtractMPvalues(value, range, nf, args.pbi.m_vnlosubv);
    else if (subkey == "ME_Generator")       ExtractMPvalues(value, range, nf, args.pbi.m_vmegen);
    else if (subkey == "RS_ME_Generator")    ExtractMPvalues(value, range, nf, args.pbi.m_vrsmegen);
    else if (subkey == "Loop_Generator")     ExtractMPvalues(value, range, nf, args.pbi.m_vloopgen);
    else if (subkey == "Integrator")         ExtractMPvalues(value, range, nf, args.pbi.m_vint);
    else if (subkey == "RS_Integrator")      ExtractMPvalues(value, range, nf, args.pbi.m_vrsint);
    else if (subkey == "PSI_ItMin")          ExtractMPvalues(value, range, nf, args.pbi.m_vitmin);
    else if (subkey == "RS_PSI_ItMin")       ExtractMPvalues(value, range, nf, args.pbi.m_vrsitmin);
  }
}

void Matrix_Element_Handler::BuildDecays
(Subprocess_Info &ACFS,const std::vector<std::string> &dectags)
{
  for (size_t i(0);i<dectags.size();++i) {
    std::string dec(dectags[i]);
    int osf=0;
    if (dec[0]=='Z') {
      dec=dec.substr(1);
      osf=1;
    }
    size_t pos(dec.find("->"));
    if (pos==std::string::npos) continue;
    Subprocess_Info ACDIS, ACDFS;
    std::string ini(dec.substr(0,pos));
    std::string fin(dec.substr(pos+2));
    ExtractFlavours(ACDIS,ini);
    ExtractFlavours(ACDFS,fin);
    if (ACDIS.m_ps.empty() || ACDFS.m_ps.empty())
      THROW(fatal_error,"Wrong decay specification");
    Subprocess_Info &CIS(ACDIS.m_ps.front());
    size_t oldsize(ACFS.m_ps.size()), cdsize(ACDFS.m_ps.size());
    ACFS.m_ps.resize(oldsize*cdsize);
    for (size_t cfss(1);cfss<cdsize;++cfss)
      for (size_t acfsi(0);acfsi<oldsize;++acfsi)
	ACFS.m_ps[cfss*oldsize+acfsi]=
	  ACFS.m_ps[(cfss-1)*oldsize+acfsi];
    for (size_t acfsi(0);acfsi<oldsize;++acfsi) {
      for (size_t cfss(0);cfss<cdsize;++cfss) {
	Subprocess_Info &CFS(ACDFS.m_ps[cfss]);
	msg_Debugging()<<METHOD<<"(): Init decay {\n"<<CIS<<CFS<<"}\n";
	if (CIS.NExternal()!=1 || CFS.NExternal()<2)
	  THROW(fatal_error,"Wrong number of particles in decay");
	if (!ACFS.m_ps[cfss*oldsize+acfsi].AddDecay(CIS,CFS,osf))
	  THROW(fatal_error,"No match for decay "+dec);
      }
    }
  }
}

void Matrix_Element_Handler::LimitCouplings
(MPSV_Map &pbi,const size_t &nfs,const std::string &pnid,
 std::vector<double> &mincpl,std::vector<double> &maxcpl,const int mode)
{
  std::string ds;
  if (!GetMPvalue(pbi,nfs,pnid,ds)) return;
  while (ds.find("*")!=std::string::npos) ds.replace(ds.find("*"),1,"-1");
  std::vector<double> cpl(ToVector<double>(ds));
  if (mode&1) {
    if (cpl.size()>mincpl.size()) mincpl.resize(cpl.size(),0);
    for (size_t i(0);i<mincpl.size();++i)
      if (cpl[i]>=0 && cpl[i]>mincpl[i]) mincpl[i]=cpl[i];
  }
  if (mode&2) {
    if (cpl.size()>maxcpl.size()) maxcpl.resize(cpl.size(),99);
    for (size_t i(0);i<maxcpl.size();++i)
      if (cpl[i]>=0 && cpl[i]<maxcpl[i]) maxcpl[i]=cpl[i];
  }
}

void Matrix_Element_Handler::BuildSingleProcessList(
  Single_Process_List_Args& args)
{
  int aoqcd(0), loprocs(0);
  Subprocess_Info AIS, AFS;
  ExtractFlavours(AIS,args.ini);
  ExtractFlavours(AFS,args.fin);
  std::vector<Process_Base*> procs;
  NLOTypeStringProcessMap_Map *pmap(NULL);
  for (size_t fss(0);fss<AFS.m_ps.size();++fss) {
    Subprocess_Info ACFS;
    ACFS.m_ps.push_back(AFS.m_ps[fss]);
    BuildDecays(ACFS,args.dectags);
    for (size_t afsi(0);afsi<ACFS.m_ps.size();++afsi) {
      msg_Debugging()<<METHOD<<"(): Check N_max ("
		     <<fss<<"): {\n"<<ACFS.m_ps[afsi]<<"}\n";
      args.pi.m_fi.GetNMax(ACFS.m_ps[afsi]);
    }
  }
  for (size_t iss(0);iss<AIS.m_ps.size();++iss) {
    Subprocess_Info &IS(AIS.m_ps[iss]);
    for (size_t fss(0);fss<AFS.m_ps.size();++fss) {
      Subprocess_Info &FS(AFS.m_ps[fss]);
      msg_Debugging()<<METHOD<<"(): Init core ("<<iss<<","
		     <<fss<<"): {\n"<<IS<<FS<<"}\n";
      std::vector<Flavour> flavs;
      IS.GetExternal(flavs);
      if (flavs.size()>1) {
        if (!p_isr->CheckConsistency(&flavs.front())) {
          msg_Error()<<METHOD<<"(): Error in initialising ISR ("
                     <<p_isr->Flav(0)<<" -> "<<flavs[0]<<") x ("
                     <<p_isr->Flav(1)<<" -> "<<flavs[1]
                     <<"). Ignoring process."<<std::endl;
          continue;
        }
      }
      Subprocess_Info ACFS;
      ACFS.m_ps.push_back(FS);
      BuildDecays(ACFS,args.dectags);
      for (size_t afsi(0);afsi<ACFS.m_ps.size();++afsi) {
	Subprocess_Info &CFS(ACFS.m_ps[afsi]);
	msg_Debugging()<<METHOD<<"(): Init process ("<<iss<<","
		       <<fss<<"): {\n"<<IS<<CFS<<"}\n";
	std::vector<Flavour> flavs;
	IS.GetExternal(flavs);
	CFS.GetExternal(flavs);
	size_t nis(IS.NExternal()), nfs(CFS.NExternal());
	double inisum=0.0, finsum=0.0, dd(0.0);
	for (size_t i(0);i<nis;++i) inisum+=flavs[i].Mass();
	for (size_t i(0);i<nfs;++i) finsum+=flavs[i+nis].Mass();
	if (inisum>rpa->gen.Ecms() || finsum>rpa->gen.Ecms()) continue;
	std::string pnid(CFS.MultiplicityTag()), ds;
	int di;
	Process_Info cpi(args.pi);
	cpi.m_ii=IS;
	cpi.m_fi=CFS;
	cpi.m_fi.m_nlotype=args.pi.m_fi.m_nlotype;
	cpi.m_fi.m_nlocpl=args.pi.m_fi.m_nlocpl;
	cpi.m_fi.SetNMax(args.pi.m_fi);
	LimitCouplings(args.pbi.m_vmincpl,nfs,pnid,cpi.m_mincpl,cpi.m_maxcpl,1);
	LimitCouplings(args.pbi.m_vmaxcpl,nfs,pnid,cpi.m_mincpl,cpi.m_maxcpl,2);
	LimitCouplings(args.pbi.m_vcpl,nfs,pnid,cpi.m_mincpl,cpi.m_maxcpl,3);
	LimitCouplings(args.pbi.m_vminacpl,nfs,pnid,cpi.m_minacpl,cpi.m_maxacpl,1);
	LimitCouplings(args.pbi.m_vmaxacpl,nfs,pnid,cpi.m_minacpl,cpi.m_maxacpl,2);
	LimitCouplings(args.pbi.m_vacpl,nfs,pnid,cpi.m_minacpl,cpi.m_maxacpl,3);
	// automatically increase QCD coupling for QCD multijet merging
	if (cpi.m_ckkw&1) {
	  cpi.m_mincpl[0]+=aoqcd;
	  cpi.m_maxcpl[0]+=aoqcd;
	  cpi.m_minacpl[0]+=aoqcd;
	  cpi.m_maxacpl[0]+=aoqcd;
	  ++aoqcd;
	}

	// test whether cpls are halfinteger, fill in open spots for same size
	size_t minsize(Min(cpi.m_mincpl.size(),cpi.m_maxcpl.size()));
	double intpart,fracpart;
	for (size_t i(0);i<cpi.m_mincpl.size();++i) {
	  fracpart=std::modf(2.*cpi.m_mincpl[i],&intpart);
	  if (fracpart!=0.)
	    THROW(fatal_error,"Order/Min_Order contains non-halfinteger entry "
			      +std::string("at position ")
			      +ToString(i)+": "+ToString(intpart)+"."
			      +ToString(std::abs(fracpart)).substr(1)+". Abort.");
	}
	for (size_t i(0);i<cpi.m_maxcpl.size();++i) {
	  fracpart=std::modf(2.*cpi.m_maxcpl[i],&intpart);
	  if (fracpart!=0.)
	    THROW(fatal_error,"Order/Max_Order contains non-halfinteger entry "
			      +std::string("at position ")
			      +ToString(i)+": "+ToString(intpart)+"."
			      +ToString(fracpart)+". Abort.");
	}
	for (size_t i(0);i<minsize;++i) {
	  if (cpi.m_mincpl[i]>cpi.m_maxcpl[i]) {
	    msg_Error()<<METHOD<<"(): Invalid coupling orders: "
		       <<cpi.m_mincpl<<" .. "<<cpi.m_maxcpl<<"\n";
	    THROW(inconsistent_option,"Please correct coupling orders");
	  }
	}
	if (GetMPvalue(args.pbi.m_vscale,nfs,pnid,ds)) cpi.m_scale=ds;
	if (GetMPvalue(args.pbi.m_vcoupl,nfs,pnid,ds)) cpi.m_coupling=ds;
	if (GetMPvalue(args.pbi.m_vkfac,nfs,pnid,ds)) cpi.m_kfactor=ds;
        cpi.m_selectors = args.pbi.m_selectors;
	if (GetMPvalue(args.pbi.m_vnmaxq,nfs,pnid,di)) cpi.m_nmaxq=di;
	if (GetMPvalue(args.pbi.m_vnminq,nfs,pnid,di)) cpi.m_nminq=di;
	if (GetMPvalue(args.pbi.m_vamegicmhv,nfs,pnid,di)) cpi.m_amegicmhv=di;
	if (GetMPvalue(args.pbi.m_vntchan,nfs,pnid,di)) cpi.m_ntchan=di;
	if (GetMPvalue(args.pbi.m_vmtchan,nfs,pnid,di)) cpi.m_mtchan=di;
	if (GetMPvalue(args.pbi.m_vgpath,nfs,pnid,ds)) cpi.m_gpath=ds;
	if (GetMPvalue(args.pbi.m_vaddname,nfs,pnid,ds)) cpi.m_addname=ds;
	if (GetMPvalue(args.pbi.m_vspecial,nfs,pnid,ds)) cpi.m_special=ds;
	if (GetMPvalue(args.pbi.m_vnlomode,nfs,pnid,ds)) {
	  args.pi.m_nlomode=cpi.m_nlomode=ToType<nlo_mode::code>(ds);
	  if (cpi.m_nlomode==nlo_mode::unknown)
	    THROW(fatal_error,"Unknown NLO_Mode "+ds+" {"+pnid+"}");
	  cpi.m_fi.m_nlotype=ToType<nlo_type::code>("BVIRS");
	  if (m_nlomode==nlo_mode::none) m_nlomode=cpi.m_nlomode;
	  if (cpi.m_nlomode!=m_nlomode)
	    THROW(fatal_error,"Unable to process multiple NLO modes at the "
			      "same time");
	}
	if (GetMPvalue(args.pbi.m_vnlopart,nfs,pnid,ds)) {
	  cpi.m_fi.m_nlotype=ToType<nlo_type::code>(ds);
	}
	if (GetMPvalue(args.pbi.m_vnlocpl,nfs,pnid,ds)) {
          cpi.m_fi.m_nlocpl = ToVector<double>(ds);
	  if (cpi.m_fi.m_nlocpl.size()<2) cpi.m_fi.m_nlocpl.resize(2,0);
	}
	if (GetMPvalue(args.pbi.m_vnlosubv,nfs,pnid,ds)) cpi.m_fi.m_sv=ds;
	if (GetMPvalue(args.pbi.m_vmegen,nfs,pnid,ds)) cpi.m_megenerator=ds;
	if (GetMPvalue(args.pbi.m_vrsmegen,nfs,pnid,ds)) cpi.m_rsmegenerator=ds;
	else cpi.m_rsmegenerator=cpi.m_megenerator;
	if (GetMPvalue(args.pbi.m_vloopgen,nfs,pnid,ds)) {
	  m_gens.LoadGenerator(ds);
	  cpi.m_loopgenerator=ds;
	}
	if (GetMPvalue(args.pbi.m_vint,nfs,pnid,ds)) cpi.m_integrator=ds;
	if (GetMPvalue(args.pbi.m_vrsint,nfs,pnid,ds)) cpi.m_rsintegrator=ds;
	else cpi.m_rsintegrator=cpi.m_integrator;
	if (GetMPvalue(args.pbi.m_vitmin,nfs,pnid,di)) cpi.m_itmin=di;
	if (GetMPvalue(args.pbi.m_vrsitmin,nfs,pnid,di)) cpi.m_rsitmin=di;
	else cpi.m_rsitmin=cpi.m_itmin;
	std::vector<Process_Base*> proc=InitializeProcess(cpi,pmap);
	for (size_t i(0);i<proc.size();i++) {
	  if (proc[i]==NULL)
	    msg_Error()<<METHOD<<"(): No process for {\n"
		       <<cpi<<"\n}"<<std::endl;
	  procs.push_back(proc[i]);
	  proc[i]->Integrator()->
	    SetISRThreshold(ATOOLS::Max(inisum,finsum));
	  if (GetMPvalue(args.pbi.m_vefac,nfs,pnid,dd))
	    proc[i]->Integrator()->SetEnhanceFactor(dd);
	  if (GetMPvalue(args.pbi.m_vmaxeps,nfs,pnid,dd))
	    proc[i]->Integrator()->SetMaxEpsilon(dd);
	  else proc[i]->Integrator()->SetMaxEpsilon(1.0e-3);
	  if (GetMPvalue(args.pbi.m_vrsefac,nfs,pnid,dd))
	    proc[i]->Integrator()->SetRSEnhanceFactor(dd);
	  double maxerr(-1.0);
	  std::string eobs, efunc;
	  if (GetMPvalue(args.pbi.m_vmaxerr,nfs,pnid,dd)) maxerr=dd;
	  if (GetMPvalue(args.pbi.m_veobs,nfs,pnid,ds)) eobs=ds;
	  if (GetMPvalue(args.pbi.m_vefunc,nfs,pnid,ds)) efunc=ds;
	  proc[i]->InitPSHandler(maxerr,eobs,efunc);
	  proc[i]->SetShower(p_shower->GetShower());
	}
	if (loprocs==0) loprocs=procs.size();
      }
    }
  }
  for (size_t i(0);i<procs.size();++i) {
    Process_Info &cpi(procs[i]->Info());
    Selector_Key skey;
    if (cpi.m_selectors.GetItemsCount() == 0) {
      skey.m_settings = Settings::GetMainSettings()["SELECTORS"];
    } else {
      skey.m_settings = cpi.m_selectors;
    }
    if (args.pi.m_ckkw&1) {
      MyStrStream jfyaml;
      jfyaml << "METS: {";
      std::string ycut{ args.pbi.m_gycut };
      GetMPvalue(args.pbi.m_vycut,
                 cpi.m_fi.NExternal(),
		 cpi.m_fi.MultiplicityTag(),
                 ycut);
      jfyaml << "YCUT: \"" << ycut << "\"";
      if (i<loprocs) {
        jfyaml << ", LO: true";
	if (args.pbi.m_cutcore==true) {
          jfyaml << ", CUT: true";
	}
      }
      jfyaml << "}";
      skey.AddSelectorYAML(jfyaml.str());
    }
    procs[i]->SetSelector(skey);
    procs[i]->SetScale
      (Scale_Setter_Arguments(p_model,cpi.m_scale,cpi.m_coupling));
    procs[i]->SetKFactor
      (KFactor_Setter_Arguments(cpi.m_kfactor));
  }
}

size_t Matrix_Element_Handler::ExtractFlavours(Subprocess_Info &info,
                                               std::string buffer)
{
  info.m_ps.resize(1);
  info.m_ps.front().m_ps.clear();
  while(true) {
    while (buffer.length()>0 && 
	   (buffer[0]==' ' || buffer[0]=='\t')) buffer.erase(0,1);
    if (buffer.length()==0) break;
    size_t pos(Min(buffer.find(' '),buffer.length()));
    std::string cur(buffer.substr(0,pos));
    buffer=buffer.substr(pos);
    pos=cur.find('(');
    std::string polid, mpl, rem;
    if (pos!=std::string::npos) {
      polid=cur.substr(pos);
      rem=polid.substr(polid.find(')')+1);
      cur=cur.substr(0,pos);
      polid=polid.substr(1,polid.find(')')-1);
    }
    if (cur.length()==0 && polid.length()) {
      cur="0"+rem;
      mpl=polid;
      polid="";
    }
    pos=cur.find('[');
    std::string decid;
    if (pos!=std::string::npos) {
      decid=cur.substr(pos);
      cur=cur.substr(0,pos);
      decid=decid.substr(1,decid.find(']')-1);
    }
    int n(-1);
    pos=cur.find('{');
    if (pos!=std::string::npos) {
      std::string nid(cur.substr(pos));
      cur=cur.substr(0,pos);
      n=ToType<size_t>(nid.substr(1,nid.find('}')-1));
    }
    int kfc(ToType<int>(cur));
    Flavour cfl((kf_code)abs(kfc));
    if (kfc<0) cfl=cfl.Bar();
    if (n==-1) {
      for (size_t i(0);i<info.m_ps.size();++i)
	info.m_ps[i].m_ps.push_back(Subprocess_Info(cfl,decid,polid,mpl));
    }
    else {
      size_t oldsize(info.m_ps.size());
      info.m_ps.resize(oldsize*(n+1));
      for (int j(1);j<=n;++j) {
	for (size_t i(0);i<oldsize;++i) {
	  info.m_ps[j*oldsize+i]=info.m_ps[(j-1)*oldsize+i];
	  info.m_ps[j*oldsize+i].m_ps.push_back(Subprocess_Info(cfl,"",polid,mpl));
	}
      }
    }
  }
  return info.m_ps.back().m_ps.size();
}

namespace SHERPA {

  template <typename Type>
  void Matrix_Element_Handler::AddMPvalue
  (std::string lstr,std::string rstr,const Type &val,
   std::map<std::string,std::pair<int,Type> >& dv,
   const int nfs,const int &priority)
  {
    if (rstr.length()==0) {
      if (nfs==0 && 
	  (dv.find(lstr)==dv.end() || dv[lstr].first>priority)) {
	msg_Debugging()<<METHOD<<"(): adding '"<<val
		       <<"' {"<<lstr<<"}("<<priority<<")\n";
	dv[lstr]=std::pair<int,Type>(priority,val);
      }
      return;
    }
    size_t pos(rstr.find('-')), ltp(rstr.find('['));
    if (pos==std::string::npos || ltp<pos-1) {
      if (ltp!=std::string::npos) {
	size_t rtp(rstr.find(']',ltp));
	AddMPvalue(lstr+rstr.substr(0,rtp+1),rstr.substr(rtp+1),val,dv,
		   nfs-ToType<int>(rstr.substr(ltp+1,rtp-ltp-1)),priority);
	return;
      }
      AddMPvalue(lstr+rstr,"",val,dv,nfs-ToType<int>(rstr),priority);
      return;
    }
    std::string rlstr(rstr.substr(0,pos)), rrstr(rstr.substr(pos+1)), rmstr;
    if (pos>0 && ltp==pos-1) {
      rmstr="]";
      rrstr=rrstr.substr(1);
    }
    for (int i(0);i<=nfs;++i)
      AddMPvalue(lstr+rlstr+ToString(i)+rmstr,rrstr,val,dv,nfs-i,priority);
  }

  template void Matrix_Element_Handler::AddMPvalue
  (std::string lstr,std::string rstr,const double &val,
   std::map<std::string,std::pair<int,double> >& dv,
   const int nfs,const int &priority);
  template void Matrix_Element_Handler::AddMPvalue
  (std::string lstr,std::string rstr,const std::string &val,
   std::map<std::string,std::pair<int,std::string> >& dv,
   const int nfs,const int &priority);

  template <typename Type>
  bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,Type> >& dv,
   const int nfs,const std::string &pnid,Type &rv)
  {
    std::map<std::string,std::pair<int,Type> > cdv(dv);
    for (typename std::map<std::string,std::pair<int,Type> >::const_iterator 
	   dit(dv.begin());dit!=dv.end();++dit) { 
      AddMPvalue<Type>("",dit->first,dit->second.second,
		       dv,nfs,dit->second.first);
    }
    if (dv.find(pnid)!=dv.end()) {
      rv=dv[pnid].second;
      return true;
    }
    std::string nfstag(ToString(nfs));
    if (dv.find(nfstag)!=dv.end()) {
      rv=dv[nfstag].second;
      return true;
    }
    return false;
  }

  template bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,int> >& dv,
   const int nfs,const std::string &pnid,int &rv);
  template bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,double> >& dv,
   const int nfs,const std::string &pnid,double &rv);
  template bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,std::string> >& dv,
   const int nfs,const std::string &pnid,std::string &rv);

  template <typename Type>
  void Matrix_Element_Handler::ExtractMPvalues(
    std::string str,
    std::pair<size_t, size_t> multirange,
    const int &priority,
    std::map<std::string,std::pair<int,Type> >& dv) const
  {
    if (str == "") return;
    const auto value = ToType<Type>(str);
    if (multirange.second == std::numeric_limits<size_t>::max()) {
      dv["-"] = std::pair<int, Type>(priority, value);
      msg_Debugging()<<METHOD<<"(): adding '"<<value<<"'("<<dv["-"].second
		     <<") {-}("<<priority<<")\n";
      return;
    }
    for (auto m = multirange.first; m <= multirange.second; ++m) {
      dv[ToString(m)] = std::pair<int, Type>(priority, value);
      msg_Debugging()<<METHOD<<"(): adding '"<<value
                     <<"'("<<dv[ToString(m)].second
                     <<") {"<<ToString(m)<<"}("<<priority<<")\n";
    }
  }

}

std::string Matrix_Element_Handler::MakeString
(const std::vector<std::string> &in,const size_t &first) const
{
  std::string out(in.size()>first?in[first]:"");
  for (size_t i(first+1);i<in.size();++i) out+=" "+in[i];
  return out;
}

std::string Matrix_Element_Handler::MakeOrderString(Scoped_Settings&& s) const
{
  auto orderkeys = s.GetKeys();
  std::vector<std::string> ordervalues(orderkeys.size(), "-1");
  for (auto orderkey : orderkeys) {
    const auto order = s[orderkey]
      .SetDefault("-1")
      .SetReplacementList(String_Map{{"Any", "-1"}})
      .Get<std::string>();
    const auto orderidx = p_model->IndexOfOrderKey(orderkey);
    if (orderidx + 1 > ordervalues.size())
      ordervalues.resize(orderidx + 1, "-1");
    ordervalues[orderidx] = order;
  }
  // translate back into a single string to use ExtractMPvalues below
  return MakeString(ordervalues, 0);
}

double Matrix_Element_Handler::GetWeight
(const Cluster_Amplitude &ampl,
 const nlo_type::code type,const int mode) const
{
  std::string name(Process_Base::GenerateName(&ampl));
  for (int i(0);i<m_pmaps.size();++i) {
    StringProcess_Map::const_iterator pit
      (m_pmaps[i]->find(type)->second->find(name));
    if(pit==m_pmaps[i]->find(type)->second->end()) continue;
    auto ci = pit->second->Integrator()->ColorIntegrator();
    if (ci != nullptr) {
      ci->GeneratePoint();
      for (size_t j(0);j<ampl.Legs().size();++j)
	ampl.Leg(j)->SetCol(ColorID(ci->I()[j],ci->J()[j]));
      if (mode&1) ci->SetWOn(false);
      double res(pit->second->Differential(ampl,Variations_Mode::nominal_only));
      ci->SetWOn(true);
      return res;
    }
  }
  return 0.0;
}

