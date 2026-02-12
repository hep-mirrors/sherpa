#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/Single_Events/Analysis_Phase.H"
#include "SHERPA/Single_Events/Userhook_Phase.H"
#include "SHERPA/Single_Events/Output_Phase.H"
#include "SHERPA/Single_Events/EvtReadin_Phase.H"
#include "SHERPA/Single_Events/Signal_Processes.H"
#include "SHERPA/Single_Events/Hard_Decays.H"
#include "SHERPA/Single_Events/Minimum_Bias.H"
#include "SHERPA/Single_Events/Multiple_Interactions.H"
#include "SHERPA/Single_Events/Jet_Evolution.H"
#include "SHERPA/Single_Events/Signal_Process_FS_QED_Correction.H"
#include "SHERPA/Single_Events/Beam_Remnants.H"
#include "SHERPA/Single_Events/Hadronization.H"
#include "SHERPA/Single_Events/Hadron_Decays.H"
#include "SHERPA/PerturbativePhysics/Hard_Decay_Handler.H"
#include "SHERPA/Tools/HepMC3_Interface.H"
#include "PHASIC++/Decays/Decay_Channel.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/PDF_Base.H"
#include <cstring>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Sherpa::Sherpa(int argc, char* argv[]) :
  p_inithandler(nullptr),
  p_eventhandler(nullptr)
#ifdef USING__HEPMC3
  , p_hepmc3(nullptr)
#endif
{
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  // rpa should be constructed before initializing the main settings, since the
  // latter might throw an exception and rpa would be involved in terminating
  // the program then; however, do not call its Init method yet, because this
  // in turn needs the Settings to be initialized
  ATOOLS::rpa = new Run_Parameter();
  Settings::InitializeMainSettings(argc, argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();
  m_trials = 0;
  m_debuginterval = 0;
  m_debugstep = -1;
  m_displayinterval = 100;
  m_evt_starttime = -1.0;
  exh->AddTerminatorObject(this);
}

Sherpa::~Sherpa()
{
  if (msg_LevelIsInfo()) {
    Return_Value::PrintStatistics(msg->Out());
    if (p_inithandler->GetVariations()) {
      p_inithandler->GetVariations()->PrintStatistics(msg->Out());
    }
    Blob_List::PrintMomFailStatistics(msg->Out());
    msg->PrintRates();
    PHASIC::Decay_Channel::PrintMaxKinFailStatistics(msg->Out());
  }
  if (p_eventhandler) { delete p_eventhandler; p_eventhandler = nullptr; }
  if (p_inithandler)  { delete p_inithandler;  p_inithandler  = nullptr; }
#ifdef USING__HEPMC3
  if (p_hepmc3)       { delete p_hepmc3;       p_hepmc3       = NULL; }
#endif
  Settings& s = Settings::GetMainSettings();
  if (s["CHECK_SETTINGS"].SetDefault(true).Get<bool>())
    Settings::FinalizeMainSettings();
  rpa->gen.WriteCitationInfo();
  exh->RemoveTerminatorObject(this);
  delete ATOOLS::s_loader;
  delete PDF::pdfdefs;
  delete ATOOLS::rpa;
  delete ATOOLS::ran;
#ifdef USING__MPI
  mpi->Barrier();
#endif
  delete ATOOLS::msg;
  delete ATOOLS::exh;
  delete ATOOLS::mpi;
  ATOOLS::ClearParticles();
}

bool Sherpa::InitializeTheRun()
{
  Settings& s = Settings::GetMainSettings();
  p_inithandler = new Initialization_Handler();
  RegisterDefaults();

  mpi->PrintRankInfo();

  DrawLogo(s["PRINT_VERSION_INFO"].Get<bool>());
  int initonly=s["INIT_ONLY"].Get<int>();
  if (initonly) rpa->gen.SetNumberOfEvents(0);
  if (p_inithandler->InitializeTheFramework()) {
    if (initonly==1) THROW(normal_exit,"Initialization complete.");
    if (initonly==2) return true;
    if (!p_inithandler->CalculateTheHardProcesses()) return false;
    m_showtrials=s["SHOW_NTRIALS"].Get<bool>();

    // read in from status path
    bool res(true);
    std::string statuspath(s["STATUS_PATH"].Get<std::string>());
    if (statuspath != "") {
      res=exh->ReadInStatus(statuspath);
    }

    m_debuginterval = s["DEBUG_INTERVAL"].Get<long int>();
    m_debugstep     = s["DEBUG_STEP"].Get<long int>();

    m_displayinterval=s["EVENT_DISPLAY_INTERVAL"].Get<int>();
    m_evt_output = s["EVT_OUTPUT"].Get<int>();
    m_evt_output_start = s["EVT_OUTPUT_START"].Get<int>();

    return res;
  }
  msg_Error()<<"Error in Sherpa::InitializeRun()"<<endl
	     <<"   Did not manage to initialize the framework."<<endl
	     <<"   Try to run nevertheless ... ."<<endl;

  return 0;
}

void Sherpa::RegisterDefaults()
{
  Settings& s = Settings::GetMainSettings();
  s["PRINT_VERSION_INFO"].SetDefault(false);
  s["INIT_ONLY"].SetDefault(0);
  s["SHOW_NTRIALS"].SetDefault(false);
  s["DEBUG_INTERVAL"].SetDefault(0);
  s["DEBUG_STEP"].SetDefault(-1);
  s["EVENT_DISPLAY_INTERVAL"].SetDefault(100);
  s["EVT_OUTPUT"].SetDefault(msg->Level());
  s["MSG_LIMIT"].SetDefault(20);
  msg->SetLimit(s["MSG_LIMIT"].Get<int>());

  const int evtoutput{ s["EVT_OUTPUT"].Get<int>() };
  s["EVT_OUTPUT_START"].SetDefault(evtoutput != msg->Level() ? 1 : 0);
}

bool Sherpa::InitializeTheEventHandler()
{
  eventtype::code mode = p_inithandler->Mode();
  p_eventhandler  = new Event_Handler();
  p_eventhandler->SetVariations(p_inithandler->GetVariations());
  Analysis_Vector *anas(p_inithandler->GetAnalyses());
  for (Analysis_Vector::iterator it=anas->begin(); it!=anas->end(); ++it) {
    (*it)->SetEventHandler(p_eventhandler);
  }

  if (mode==eventtype::EventReader) {
    p_eventhandler->AddEventPhase(new EvtReadin_Phase(p_inithandler->GetEventReader()));
    p_eventhandler->AddEventPhase(new Hard_Decays(p_inithandler->GetHardDecayHandler()));
    p_eventhandler->AddEventPhase(new Beam_Remnants(p_inithandler->GetBeamRemnantHandler()));
  }
  else {
    p_eventhandler->AddEventPhase(new Signal_Processes(p_inithandler->GetMatrixElementHandler()));
    p_eventhandler->AddEventPhase(new Minimum_Bias(p_inithandler->GetSoftCollisionHandlers()));
    p_eventhandler->AddEventPhase(new Hard_Decays(p_inithandler->GetHardDecayHandler()));
    p_eventhandler->AddEventPhase(new Jet_Evolution(p_inithandler->GetMatrixElementHandler(),
                                                    p_inithandler->GetHardDecayHandler(),
						    p_inithandler->GetHDHandler(),
						    p_inithandler->GetMIHandlers(),
						    p_inithandler->GetSoftCollisionHandlers(),
						    p_inithandler->GetShowerHandlers(),
						    p_inithandler->GetRemnantHandlers()));
    p_eventhandler->AddEventPhase(new Signal_Process_FS_QED_Correction(
						    p_inithandler->GetMatrixElementHandler(),
						    p_inithandler->GetSoftPhotonHandler()));
    p_eventhandler->AddEventPhase(new Multiple_Interactions(p_inithandler->GetMIHandlers()));
    p_eventhandler->AddEventPhase(new Beam_Remnants(p_inithandler->GetBeamRemnantHandler()));
    p_eventhandler->AddEventPhase(new Hadronization(p_inithandler->GetColourReconnectionHandler(),
						    p_inithandler->GetFragmentation()));
    p_eventhandler->AddEventPhase(new Hadron_Decays(p_inithandler->GetHDHandler()));
  }
  p_eventhandler->AddEventPhase(new Userhook_Phase(this));
  if (!anas->empty()) p_eventhandler->AddEventPhase(new Analysis_Phase(anas));
  if (!p_inithandler->GetOutputs()->empty())
    p_eventhandler->AddEventPhase(new Output_Phase(p_inithandler->GetOutputs(), p_eventhandler));
  p_eventhandler->SetFilter(p_inithandler->GetFilter());
  p_eventhandler->PrintGenericEventStructure();

  ran->EraseLastIncrementedSeed();

  return 1;
}


bool Sherpa::GenerateOneEvent(bool reset)
{
  if (m_evt_output_start>0 &&
      m_evt_output_start==rpa->gen.NumberOfGeneratedEvents()+1) {
    msg->SetLevel(m_evt_output);
  }

  if(m_debuginterval>0 &&
     rpa->gen.NumberOfGeneratedEvents()%m_debuginterval==0 &&
     (p_inithandler->GetMatrixElementHandler()->SeedMode()!=3 ||
      rpa->gen.NumberOfGeneratedEvents()==0)) {
      std::string fname=ToString(rpa->gen.NumberOfGeneratedEvents())+".dat";
      ran->WriteOutStatus(("random."+fname).c_str());
  }
  if (m_debugstep>=0) {
    if (p_inithandler->GetMatrixElementHandler()->SeedMode()!=3)
      ran->ReadInStatus(("random."+ToString(m_debugstep)+".dat").c_str());
    else {
      ran->ReadInStatus("random.0.dat");
      ran->FastForward(m_debugstep);
    }
  }

  if (m_evt_starttime<0.0) m_evt_starttime=rpa->gen.Timer().RealTime();

  if (reset) p_eventhandler->Reset();
  rpa->gen.SetIsGen(false);
  if (p_eventhandler->GenerateEvent(p_inithandler->Mode())) {
    if(m_debuginterval>0 && rpa->gen.NumberOfGeneratedEvents()%m_debuginterval==0){
      std::string fname=ToString(rpa->gen.NumberOfGeneratedEvents())+".dat";
      std::ofstream eventout(("refevent."+fname).c_str());
      eventout<<"# trial "<<rpa->gen.NumberOfTrials()-1<<std::endl;
      eventout<<*p_eventhandler->GetBlobs()<<std::endl;
      eventout.close();
    }
    if (m_debugstep>=0) {
      std::ofstream event(("event."+ToString(m_debugstep)+".dat").c_str());
      event<<*p_eventhandler->GetBlobs()<<std::endl;
      event.close();
      THROW(normal_exit,"Debug event written.");
    }
    rpa->gen.SetNumberOfGeneratedEvents(rpa->gen.NumberOfGeneratedEvents()+1);
    Blob_List *blobs(p_eventhandler->GetBlobs());

    /// Increase m_trials --- based on signal blob["Trials"] if existent
    if (blobs->FindFirst(btp::Signal_Process) == nullptr) {
      m_trials+=1;
      msg_Debugging()<<"  No Signal_Process Blob found, increasing m_trials by 1\n";
    }
    else {
      m_trials+=(*blobs->FindFirst(btp::Signal_Process))["Trials"]->Get<double>();
      std::string sub_name = rpa->gen.GetIsGenSubName();
      std::chrono::high_resolution_clock::time_point end3 = std::chrono::high_resolution_clock::now();
      double finetime = std::chrono::duration_cast<std::chrono::nanoseconds>(end3-rpa->gen.GetIsGenTime()).count()/1000000000.;
      rpa->gen.SetNumberMap("n_kept_"+sub_name, rpa->gen.NumberMap("n_kept_"+sub_name)+1);
      rpa->gen.SetTimeMap("sum_overhead_after_kept_" +sub_name, rpa->gen.TimeMap("sum_overhead_after_kept_" +sub_name)+finetime);
      rpa->gen.SetTimeMap("sum2_overhead_after_kept_"+sub_name, rpa->gen.TimeMap("sum2_overhead_after_kept_"+sub_name)+finetime*finetime);
      rpa->gen.SetNumberMap("n_overhead_after_kept_" +sub_name, rpa->gen.NumberMap("n_overhead_after_kept_" +sub_name)+1);
    }

    if (msg_LevelIsEvents()) {
      if (!blobs->empty()) {
	msg_Out()<<"  -------------------------------------------------\n";
	for (Blob_List::iterator blit=blobs->begin();
	     blit!=blobs->end();++blit)
	  msg_Out()<<*(*blit)<<std::endl;
	msg_Out()<<"  -------------------------------------------------\n";
      }
      else msg_Out()<<"  ******** Empty event ********  "<<std::endl;
    }

    int i=rpa->gen.NumberOfGeneratedEvents();
    int nevt=rpa->gen.NumberOfEvents();
    msg_Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
    int exp;
    for (exp=5; i/int(pow(10,exp))==0; --exp) {}
    if (((rpa->gen.BatchMode()&4 && i%m_displayinterval==0) ||
	 (!(rpa->gen.BatchMode()&4) && i%int(pow(10,exp))==0)) &&
	i<rpa->gen.NumberOfEvents()) {
      double diff=rpa->gen.Timer().RealTime()-m_evt_starttime;
      msg_Info()<<"  Event "<<i;
      if (m_showtrials)
        msg_Info()<<"("+ToString(m_trials)+")";
      msg_Info()<<" ( ";
      if (rpa->gen.BatchMode()&16) {
        msg_Info()<<diff<<"s elapsed / "
                  <<((nevt-i)/(double)i*diff)<<"s";
      } else {
        msg_Info()<<FormatTime(size_t(diff))<<" elapsed / "
                  <<FormatTime(size_t((nevt-i)/(double)i*diff));
      }
      msg_Info()<<" left ) -> ETA: "<<rpa->gen.Timer().
        StrFTime("%a %b %d %H:%M",time_t((nevt-i)/(double)i*diff))<<"  ";
      p_eventhandler->PerformMemoryMonitoring();
      const Uncertain<double> xs = p_eventhandler->TotalNominalXSMPI();
      if (!(rpa->gen.BatchMode()&2)) msg_Info()<<"\n  ";
      msg_Info() << "XS = " << xs.value << " pb +- ( " << xs.error
                 << " pb = " << xs.PercentError() << " % )  ";
      if (rpa->gen.BatchMode()&8)
        msg_Info()<<"  Process was "<<p_eventhandler->CurrentProcess()<<"  ";
      if (!(rpa->gen.BatchMode()&2))
	msg_Info()<<mm(1,mm::up);
      if (rpa->gen.BatchMode()&2) { msg_Info()<<std::endl; }
      else { msg_Info()<<bm::cr<<std::flush; }
    }
    return 1;
  }
  return 0;
}

#ifdef USING__HEPMC3
void Sherpa::FillHepMCEvent(HepMC3::GenEvent& event)
{
  if (p_hepmc3==NULL) p_hepmc3 = new SHERPA::HepMC3_Interface();
  ATOOLS::Blob_List* blobs=GetEventHandler()->GetBlobs();
  p_hepmc3->Sherpa2HepMC(blobs, event);
  p_hepmc3->AddCrossSection(event, p_eventhandler->TotalXS(),
                            p_eventhandler->TotalErr());
}
#endif

Uncertain<double> Sherpa::TotalNominalXS() const
{
  return p_eventhandler->TotalNominalXS();
}

std::string Sherpa::PDFInfo()
{
  std::string pdf="Unknown";
  PDF::ISR_Handler* isr=GetInitHandler()->GetISRHandler(PDF::isr::hard_process);
  if (isr) {
    if (isr->PDF(0)) {
      pdf=isr->PDF(0)->Type();
      if (isr->PDF(1) && isr->PDF(1)->Type()!=pdf) {
        pdf="Unknown";
      }
    }
  }
  return pdf;
}

void Sherpa::PrepareTerminate()
{
  SummarizeRun();
  exh->RemoveTerminatorObject(this);
}

bool Sherpa::SummarizeRun()
{
  if (p_eventhandler) {
    msg_Info()<<"  Event "<<rpa->gen.NumberOfGeneratedEvents()<<" ( "
              <<size_t(rpa->gen.Timer().RealTime()-m_evt_starttime)
              <<" s total ) = "
              << rpa->gen.NumberOfGeneratedEvents()*3600*24/
                 ((size_t) rpa->gen.Timer().RealTime()-m_evt_starttime)
              <<" evts/day ";
    //calculate sum_map etc.
    std::map<std::string, double>  time_map = rpa->gen.TimeMapAll();
    std::map<std::string, int>     number_map = rpa->gen.NumberMapAll();
    std::vector<std::string> num_types = {"total", "trial", "PS", "ME", "gen", "overw", "maxoverw", "kept"};//
    std::map<std::string, double>  sum_map;
    std::map<std::string, double>  sum_mult_map;//defined as double to simplify following divisions
    for (const std::string& num_type : num_types) {
      sum_map[num_type] = 0;
    }
    for (auto const& [key, val] : time_map) {
      if (key.rfind("sum_PS_", 0) == 0) {
        std::string sub_name = key.substr(7);
        for (const std::string& num_type : num_types) {
          sum_map[num_type] += number_map["n_"+num_type+"_"+sub_name];
          std::string mult=sub_name.substr(0,4);
          if (sum_mult_map.find(mult+"_"+num_type) == sum_mult_map.end()) {
            sum_mult_map[mult+"_"+num_type] = 0;
          }
          sum_mult_map[mult+"_"+num_type] += number_map["n_"+num_type+"_"+sub_name];
        }
      }
    }
    
    //calculate sudakov_efficiency per process (needs sum_map etc.)
    double sepsum = 0;
    std::map<std::string, double> sudakov_efficiency;
    //uncorrelated between subprocesses
    std::map<std::string, double> sudakov_efficiency_up;
    std::map<std::string, double> sudakov_efficiency_down;
    //correlated between subprocesses
    std::map<std::string, double> sudakov_efficiency_up_corr;
    std::map<std::string, double> sudakov_efficiency_down_corr;
    for (auto const& [key, val] : time_map) {
      if (key.rfind("sum_PS_", 0) == 0) {
	std::string sub_name = key.substr(7);
	if (number_map["n_kept_"+sub_name]>0) {//take sub specific
	  sudakov_efficiency[sub_name] = number_map["n_kept_"+sub_name]*1.0/number_map["n_gen_"+sub_name];
	  double d_sud = pow(number_map["n_kept_"+sub_name]*max(1, number_map["n_gen_"+sub_name]-number_map["n_kept_"+sub_name])*1.0/pow(number_map["n_gen_"+sub_name],3),0.5);
	  sudakov_efficiency_up[sub_name] = sudakov_efficiency[sub_name]+d_sud;
	  sudakov_efficiency_down[sub_name] = sudakov_efficiency[sub_name]-d_sud;
	  sudakov_efficiency_up_corr[sub_name] = sudakov_efficiency[sub_name];
	  sudakov_efficiency_down_corr[sub_name] = sudakov_efficiency[sub_name];
	} else {
	  std::string mult=sub_name.substr(0,4);
	  if (sum_mult_map.find(mult+"_kept") != sum_mult_map.end() and sum_mult_map[mult+"_kept"]>0) {
	    //assume average of same multi: e.g. 2_6
	    sudakov_efficiency[sub_name] = sum_mult_map[mult+"_kept"]/sum_mult_map[mult+"_gen"];
	  } else {
	    //assume average of whole sample
	    sudakov_efficiency[sub_name] = sum_map["kept"]/sum_map["gen"];
	  }
	  sudakov_efficiency_up[sub_name] = sudakov_efficiency[sub_name];
	  sudakov_efficiency_down[sub_name] = sudakov_efficiency[sub_name];
	  sudakov_efficiency_up_corr[sub_name] = 1.0;
	  sudakov_efficiency_down_corr[sub_name] = 0.0;
	}
	//msg_Info() << sub_name << ": sudakov_efficiency[sub_name]: " << sudakov_efficiency[sub_name] << std::endl;
	double this_sepsum = time_map["sum_overhead_after_"+sub_name]+time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_total_"+sub_name];
	sepsum += this_sepsum;
      }
    }
    //calculate chosen effevperev (needs sudakov_efficiency)
    std::map<std::string, double> chosen_alpha_map = rpa->gen.AlphaMap();
    std::map<std::string, double> chosen_efficiency_map = rpa->gen.EfficiencyMap();
    std::map<std::string, double> xsec_map = rpa->gen.XsecMap();
    double sum_p_unw = 0;
    double sum_p_eff = 0;
    double sum_p_eff_sign = 0;
    for (auto const& [key, val] : chosen_alpha_map) {
      std::string sub_name = key;
      double curr_xsec = dabs(xsec_map[sub_name])/chosen_efficiency_map[sub_name]/chosen_alpha_map[sub_name];
      sum_p_unw += chosen_efficiency_map[sub_name]*sudakov_efficiency[sub_name]*curr_xsec;
      sum_p_eff += chosen_efficiency_map[sub_name]*sudakov_efficiency[sub_name]*chosen_alpha_map[sub_name]*curr_xsec;
      sum_p_eff_sign += xsec_map[sub_name]*sudakov_efficiency[sub_name];
    }
    double chosen_effevperev = sum_p_eff/sum_p_unw*pow(sum_p_eff_sign/sum_p_eff,2);
    int generation_mode=ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE"));
    if (generation_mode!=0) {
      msg_Info()<<"with "<< chosen_effevperev << " Neff/evt           "<<std::endl;
    } else {
      //for weighted the chosen efficiency is not correctly set yet - one can see it in epsilon_scan. If weighted stays implemented with its own selectionWeight, then one could add it here later
      msg_Info()<<"                      "<<std::endl;
    }
    p_eventhandler->Finish();

    Settings& s = Settings::GetMainSettings();
    double timing_statistics_large_weight_fraction=s["TIMING_STATISTICS_LARGE_WEIGHT_FRACTION"].SetDefault(0.001).Get<double>();
    double timing_statistics_det_sim=s["TIMING_STATISTICS_DET_SIM_IN_S"].SetDefault(0.0).Get<double>();
    int timing_statistics=s["TIMING_STATISTICS"].SetDefault(0).Get<int>();
    if (not timing_statistics) return true;
    double m_ovwth = s["OVERWEIGHT_THRESHOLD"].SetDefault(1e12).Get<double>();

    //loop over subprocesses
    std::vector<std::string> sum_types = {"total", "PS", "ME", "overhead", "overhead_after_kept", "overhead_after"};//
    std::map<std::string, double> total;
    for (const std::string& sum_type : sum_types) {
      total[sum_type] = 0;
    }
    for (auto const& [key, val] : time_map) {
      if (key.rfind("sum_PS_", 0) == 0) {
        std::string sub_name = key.substr(7);
        if (timing_statistics>4) {
	  msg_Info() << std::endl;
	  for (const std::string& num_type : num_types) {
	    msg_Info()<<sub_name<<" : "<< "n_"<<num_type<<": "<<number_map["n_"+num_type+"_"+sub_name]<<std::endl;
	  }
	  for (const std::string& sum_type : sum_types) {
	    if (sum_type=="overhead") continue;
	    printTime(sub_name, sum_type, time_map, number_map);
	  }
	  printEffi(sub_name, "Cut", "ME", "trial", number_map);
	  printEffi(sub_name, "Unweighting", "gen", "ME", number_map);
	  printEffi(sub_name, "Total", "gen", "trial", number_map);
	  double w = time_map["sum_PS_"+sub_name]+time_map["sum_ME_"+sub_name];
	  msg_Info()<< sub_name << " : "<<"Summed time: "<<w<<" s"<<std::endl;
        }
        total["PS"] += time_map["sum_PS_"+sub_name];
        total["ME"] += time_map["sum_ME_"+sub_name];
        total["overhead"] += time_map["sum_total_"+sub_name]-(time_map["sum_PS_"+sub_name]+time_map["sum_ME_"+sub_name]);

        double w = time_map["sum_total_"+sub_name];
        if (timing_statistics>4) {
          msg_Info()<< sub_name << " : "<<"Total time: "<<w<<" s"<<std::endl;
        }
        total["total"] += w;
        w = time_map["sum_overhead_after_kept_"+sub_name];
        if (timing_statistics>4) {
          msg_Info()<< sub_name << " : "<<"overhead_after_kept time: "<<w<<" s"<<std::endl;
        }
        total["overhead_after_kept"] += w;
        w = time_map["sum_overhead_after_"+sub_name];
        if (timing_statistics>4) {
          msg_Info()<< sub_name << " : "<<"overhead_after time: "<<w<<" s"<<std::endl;
        }
        total["overhead_after"] += w;
      }
    }
    if (timing_statistics>4) {
      for (const std::string& num_type : num_types) {
        msg_Info()<<"sum_map['n_"<<num_type<<"']: "<<sum_map[num_type]<<std::endl;
      }
      for ( const auto &myPair : sum_mult_map) {
        msg_Info()<<"sum_mult_map['"<<myPair.first<<"']: "<<sum_mult_map[myPair.first]<<std::endl;
      }
      for (const std::string& sum_type : sum_types) {
        msg_Info()<<"total['"<<sum_type<<"']: "<<total[sum_type]<<std::endl;
      }
    }
    //following for chosen epsilon max
    //Total time and splitting into unweighting steps
    total["overhead_during"] = total["overhead"];
    total["overhead_after"] = total["overhead_after_kept"]+total["overhead_after"];
    if (timing_statistics>4) {
      msg_Info()<<std::endl<<"Total time: "<< total["total"]+total["overhead_after"] <<" s"<<std::endl;
      msg_Info()<<" (It sums up to 100%. 'overhead' is split into 'during' and 'after' unweighting.)"<<std::endl;
      for (const std::string& sum_type : {"PS", "ME", "overhead_during", "overhead_after"}) {
        msg_Info()<<" "<< std::right<<std::setw(2) << round(total[sum_type]/(total["total"]+total["overhead_after"])*100)<<"% in '"<<sum_type<<"'"<<std::endl;
      }
    }
    if (timing_statistics>4) msg_Info() << "Total separate sum: " << sepsum << "s resulting in: " << 60*60*24/sepsum*sum_map["kept"] << "evts/day" << std::endl;

    //epsilon_max scan manual definition
    std::vector<double> epsilon_values = rpa->gen.EpsilonValues();
    std::vector<double> fraction_values = rpa->gen.EpsilonValues();
    std::map<std::string, std::vector<double>> alpha_manual_map = rpa->gen.AlphaManualMapAll();
    std::map<std::string, std::vector<double>> wmax_manual_map = rpa->gen.WmaxManualMapAll();
    std::map<std::string, std::vector<double>> efficiency_manual_map = rpa->gen.EfficiencyManualMapAll();
    std::vector<double> mean_manual_alpha(epsilon_values.size()+2, -1);
    std::vector<double> mean_manual_events(epsilon_values.size()+2, -1);
    std::vector<double> mean_manual_eff_events(epsilon_values.size()+2, -1);
    std::vector<double> mean_manual_events_up(epsilon_values.size()+1, -1);
    std::vector<double> mean_manual_events_down(epsilon_values.size()+1, -1);
    std::vector<std::map<std::string, double>> mean_manual_events_time(epsilon_values.size()+1);
    int optimal_manual_i = 0;
    double optimal_manual_sum_t_trial = 0;
    double plain_xsec_sum = 0;
    for(int i=0; i < epsilon_values.size()+1; i++){
      double sum_alpha = 0;
      double sum_t_trial = 0;
      double sum_p_unw = 0;
      double sum_p_eff = 0;
      double sum_p_eff_sign = 0;
      double sum_p_unw_up = 0;
      double sum_p_unw_down = 0;
      double sum_p_unw_up_corr = 0;
      double sum_p_unw_down_corr = 0;
      double sum_t_ME = 0;
      double sum_t_PS = 0;
      double sum_t_ov_during = 0;
      double sum_t_ov_after = 0;
      //double current_epsilon = epsilon_values[i];
      for (auto const& [key, val] : alpha_manual_map) {
	std::string sub_name = key;
	if (wmax_manual_map[sub_name][i]==-2) continue; //this means that whisto is empty
	//std::cout << sub_name << std::endl;
	double curr_xsec = dabs(xsec_map[sub_name])/efficiency_manual_map[sub_name][i]/alpha_manual_map[sub_name][i]; //need to weight with sampling probability in manual approach. Why not sudakov? - bacause happens afterwards - but still more events needed for optimal eff events? no
	if (i==0) plain_xsec_sum += dabs(xsec_map[sub_name]);
	if (alpha_manual_map[sub_name][i]==-1) {
	  msg_Info() << "WARNING: for " << sub_name << " there is no alpha value for i=" << i << " corresponding to eps=" << exp(log(10)*epsilon_values[i]) << std::endl;
	}
	double tges  = (time_map["sum_total_"+sub_name])/number_map["n_total_"+sub_name]; //in s
	if (number_map["n_total_"+sub_name]==0) {
	  tges = 0; //critical, because underestimate - need to make sure that enough events generated...unc estimate hard, because 0 is 0
	}
	double overhead_after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name])/number_map["n_gen_"+sub_name];
	if (number_map["n_gen_"+sub_name]==0) {
	  overhead_after = 0;
	}
	sum_t_trial += (tges+efficiency_manual_map[sub_name][i]*(overhead_after+sudakov_efficiency[sub_name]*timing_statistics_det_sim))*curr_xsec;
	sum_p_unw += efficiency_manual_map[sub_name][i]*sudakov_efficiency[sub_name]*curr_xsec;
	sum_p_eff += efficiency_manual_map[sub_name][i]*sudakov_efficiency[sub_name]*alpha_manual_map[sub_name][i]*curr_xsec;
	sum_p_eff_sign += xsec_map[sub_name]*sudakov_efficiency[sub_name];

	//100% correlated for very few kept events
	sum_p_unw_up_corr += efficiency_manual_map[sub_name][i]*sudakov_efficiency_up_corr[sub_name]*curr_xsec;
	sum_p_unw_down_corr += efficiency_manual_map[sub_name][i]*sudakov_efficiency_down_corr[sub_name]*curr_xsec;
	//0% correlated for a lot kept events - could even make 2 sudakov uncertainties for this ...
	sum_p_unw_up += pow(efficiency_manual_map[sub_name][i]*(sudakov_efficiency_up[sub_name]-sudakov_efficiency[sub_name])*curr_xsec,2);
	sum_p_unw_down += pow(efficiency_manual_map[sub_name][i]*(sudakov_efficiency[sub_name]-sudakov_efficiency_down[sub_name])*curr_xsec,2);

	//split by evegen step
	if (number_map["n_total_"+sub_name]>0) {
	  sum_t_ME += time_map["sum_ME_"+sub_name]/number_map["n_total_"+sub_name]*curr_xsec;
	  sum_t_PS += time_map["sum_PS_"+sub_name]/number_map["n_total_"+sub_name]*curr_xsec;
	  sum_t_ov_during += (time_map["sum_total_"+sub_name]-time_map["sum_ME_"+sub_name]-time_map["sum_PS_"+sub_name])/number_map["n_total_"+sub_name]*curr_xsec;
	}
	sum_t_ov_after += (overhead_after+sudakov_efficiency[sub_name]*timing_statistics_det_sim)*efficiency_manual_map[sub_name][i]*curr_xsec;

	sum_alpha += alpha_manual_map[sub_name][i]*efficiency_manual_map[sub_name][i]*sudakov_efficiency[sub_name]*curr_xsec; //=sum_p_eff - everything consistent :)
      }
      mean_manual_alpha[i] = sum_alpha/sum_p_unw*pow(sum_p_eff_sign/sum_p_eff,2);
      mean_manual_events[i] = 60*60*24/(sum_t_trial/sum_p_unw);
      mean_manual_eff_events[i] = mean_manual_events[i]*mean_manual_alpha[i];
      mean_manual_events_up[i] = 60*60*24/(sum_t_trial/(sum_p_unw_up_corr+pow(sum_p_unw_up,0.5)));
      mean_manual_events_down[i] = 60*60*24/(sum_t_trial/(sum_p_unw_down_corr-pow(sum_p_unw_down,0.5)));
      mean_manual_events_time[i]["ME"] = sum_t_ME;
      mean_manual_events_time[i]["PS"] = sum_t_PS;
      mean_manual_events_time[i]["ov_during"] = sum_t_ov_during;
      mean_manual_events_time[i]["ov_after"] = sum_t_ov_after;
      mean_manual_events_time[i]["sum"] = sum_t_trial;
      //"=" to update optimal_manual_sum_eff, if i=0 is optimal
      if (mean_manual_eff_events[i]>=mean_manual_eff_events[optimal_manual_i] && epsilon_values.size()>i) {
	optimal_manual_i=i;
	optimal_manual_sum_t_trial=sum_t_trial;
      }
    }
    //for weighted: take care of potentially non-optimal selection weight
    int i = epsilon_values.size()+1;
    double sum_t_trial = 0;
    sum_p_unw = 0;
    double sum_xsec = 0;
    double sum_complex = 0;
    double sum_effiselw = 0;
    for (auto const& [key, val] : alpha_manual_map) {
      std::string sub_name = key;
      if (wmax_manual_map[sub_name][i]==-2) continue; //this means that whisto is empty
      //std::cout << sub_name << std::endl;
      double selw = efficiency_manual_map[sub_name][i];//used as selw - also set like this in integrator
      double tges  = (time_map["sum_total_"+sub_name])/number_map["n_total_"+sub_name]; //in s
      if (number_map["n_total_"+sub_name]==0) {
	tges = 0; //critical, because underestimate - need to make sure that enough events generated...unc estimate hard, because 0 is 0
      }
      double overhead_after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name])/number_map["n_gen_"+sub_name];
      if (number_map["n_gen_"+sub_name]==0) {
	overhead_after = 0;
      }
      sum_t_trial += (tges+(overhead_after+sudakov_efficiency[sub_name]*timing_statistics_det_sim))*selw;
      sum_p_unw += sudakov_efficiency[sub_name]*selw;

      sum_xsec += xsec_map[sub_name]*sudakov_efficiency[sub_name];
      sum_effiselw += sudakov_efficiency[sub_name]*selw;
      //std::cout << " alpha_manual_map[sub_name][i]:" << alpha_manual_map[sub_name][i] << std::endl;
      //std::cout << " sudakov_efficiency[sub_name]:" << sudakov_efficiency[sub_name] << std::endl;
      //std::cout << " selw:" << selw << std::endl;
      //std::cout << " xsec_map[sub_name]:" << xsec_map[sub_name] << std::endl;
      //std::cout << " sum_complex+:" << pow(xsec_map[sub_name]*sudakov_efficiency[sub_name],2)/(alpha_manual_map[sub_name][i]*sudakov_efficiency[sub_name]*selw) << std::endl;
      sum_complex += pow(xsec_map[sub_name],2)*sudakov_efficiency[sub_name]/(alpha_manual_map[sub_name][i]*selw);
      //std::cout << " sum_complex:" << sum_complex << std::endl;
    }
    mean_manual_alpha[i] = pow(sum_xsec,2)/(sum_effiselw*sum_complex);//von Zettel
    mean_manual_events[i] = 60*60*24/(sum_t_trial/sum_p_unw);
    mean_manual_eff_events[i] = mean_manual_events[i]*mean_manual_alpha[i];

    //same as above, but for fraction 0.001
    //epsilon_max scan manual definition
    std::map<std::string, std::vector<double>> alpha_manual_fraction_map = rpa->gen.AlphaManualFractionMapAll();
    std::vector<double> mean_manual_fraction_alpha(epsilon_values.size()+1, -1);
    std::vector<double> mean_manual_fraction_events(epsilon_values.size()+1, -1);
    std::vector<double> mean_manual_fraction_eff_events(epsilon_values.size()+1, -1);
    int optimal_manual_fraction_i = 0;
    double optimal_manual_fraction_sum_t_trial = 0;
    for(int i=0; i < epsilon_values.size(); i++){
      double sum_complex = 0;
      double sum_t_trial = 0;
      double sum_p_unw = 0;
      double sum_p_eff_sign = 0;
      for (auto const& [key, val] : alpha_manual_fraction_map) {
	std::string sub_name = key;
	if (wmax_manual_map[sub_name][i]==-2) continue; //this means that whisto is empty
	//msg_Info() << sub_name << std::endl;
	double curr_xsec = dabs(xsec_map[sub_name])/efficiency_manual_map[sub_name][i]/alpha_manual_map[sub_name][i]; //need to weight with sampling probability in manual approach. Why not sudakov? - bacause happens afterwards - but still more events needed for optimal eff events? no
	if (alpha_manual_fraction_map[sub_name][i]==-1) {
	  msg_Info() << "WARNING: for " << sub_name << " there is no alpha value for i=" << i << " corresponding to eps=" << exp(log(10)*epsilon_values[i]) << std::endl;
	}
	double tges  = (time_map["sum_total_"+sub_name])/number_map["n_total_"+sub_name]; //in s
	if (number_map["n_total_"+sub_name]==0) {
	  tges = 0; //critical, because underestimate - need to make sure that enough events generated...unc estimate hard, because 0 is 0
	}
	double overhead_after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name])/number_map["n_gen_"+sub_name];
	if (number_map["n_gen_"+sub_name]==0) {
	  overhead_after = 0;
	}
        sum_t_trial += (tges+efficiency_manual_map[sub_name][i]*(overhead_after+sudakov_efficiency[sub_name]*timing_statistics_det_sim))*curr_xsec;
        sum_p_unw += efficiency_manual_map[sub_name][i]*sudakov_efficiency[sub_name]*curr_xsec;
        sum_p_eff_sign += xsec_map[sub_name]*sudakov_efficiency[sub_name];
	//sum_complex += pow(xsec_map[sub_name]*sudakov_efficiency[sub_name],2)/(alpha_manual_fraction_map[sub_name][i]*dabs(xsec_map[sub_name])/alpha_manual_map[sub_name][i]*sudakov_efficiency[sub_name]);
	sum_complex += dabs(xsec_map[sub_name])*sudakov_efficiency[sub_name]/(alpha_manual_fraction_map[sub_name][i]/alpha_manual_map[sub_name][i]);
      }
      mean_manual_fraction_alpha[i] = pow(sum_p_eff_sign,2)/(sum_p_unw*sum_complex);
      mean_manual_fraction_events[i] = 60*60*24/(sum_t_trial/sum_p_unw);
      mean_manual_fraction_eff_events[i] = mean_manual_fraction_events[i]*mean_manual_fraction_alpha[i];
      //"=" to update optimal_manual_sum_eff, if i=0 is optimal
      if (mean_manual_fraction_eff_events[i]>=mean_manual_fraction_eff_events[optimal_manual_fraction_i]) {
	optimal_manual_fraction_i=i;
	optimal_manual_fraction_sum_t_trial=sum_t_trial;
      }
    }





    //same as above, but for all fractions
    //fraction scan manual definition
    std::map<std::string, std::vector<std::vector<double>>> alpha_manual_fscan_map = rpa->gen.AlphaManualFScanMapAll();
    std::map<std::string, std::vector<std::vector<double>>> efficiency_manual_fscan_map = rpa->gen.EfficiencyManualFScanMapAll();
    std::vector<double> mean_manual_alpha_fscan(fraction_values.size(), -1);
    std::vector<double> mean_manual_events_fscan(fraction_values.size(), -1);
    std::vector<double> mean_manual_eff_events_fscan(fraction_values.size(), -1);
    for(int fi=0; fi < fraction_values.size(); fi++){
      double sum_alpha = 0;
      double sum_t_trial = 0;
      double sum_p_unw = 0;
      double sum_p_eff = 0;
      double sum_p_eff_sign = 0;
      for (auto const& [key, val] : alpha_manual_fscan_map) {
	std::string sub_name = key;
	//msg_Info() << "  " << sub_name << std::endl;
	double opt_alpha = -1;
	double opt_t_trial = 0;
	double opt_p_unw = 0;
	double opt_p_eff = 0;
	double opt_p_eff_sign = 0;
	for(int i=0; i < epsilon_values.size(); i++){	
	  //msg_Info() << "   " << i << std::endl;
	  double curr_xsec = dabs(xsec_map[sub_name])/efficiency_manual_fscan_map[sub_name][fi][i]/alpha_manual_fscan_map[sub_name][fi][i];
	  if (alpha_manual_fscan_map[sub_name][fi][i]==-1) {
	    msg_Info() << "WARNING: for " << sub_name << " there is no alpha value for i=" << i << " corresponding to fraction=" << exp(log(10)*fraction_values[i]) << std::endl;
	  }
	  double tges  = (time_map["sum_total_"+sub_name])/number_map["n_total_"+sub_name]; //in s
	  if (number_map["n_total_"+sub_name]==0) {
	    tges = 0; //critical, because underestimate - need to make sure that enough events generated...unc estimate hard, because 0 is 0
	  }
	  double overhead_after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name])/number_map["n_gen_"+sub_name];
	  if (number_map["n_gen_"+sub_name]==0) {
	    overhead_after = 0;
	  }
	  double this_alpha = alpha_manual_fscan_map[sub_name][fi][i]*efficiency_manual_fscan_map[sub_name][fi][i]*sudakov_efficiency[sub_name]*curr_xsec; //=sum_p_eff - everything consistent :)
	  double this_t_trial = (tges+efficiency_manual_fscan_map[sub_name][fi][i]*(overhead_after+sudakov_efficiency[sub_name]*timing_statistics_det_sim))*curr_xsec;
	  double this_p_unw = efficiency_manual_fscan_map[sub_name][fi][i]*sudakov_efficiency[sub_name]*curr_xsec;
	  double this_p_eff = efficiency_manual_fscan_map[sub_name][fi][i]*sudakov_efficiency[sub_name]*alpha_manual_fscan_map[sub_name][fi][i]*curr_xsec;
	  double this_p_eff_sign = xsec_map[sub_name]*sudakov_efficiency[sub_name];
	  if (opt_alpha==-1 or opt_t_trial/opt_p_eff/(this_t_trial/this_p_eff)>1.0001 or (opt_t_trial/opt_p_eff/(this_t_trial/this_p_eff)>0.9999 and opt_alpha/opt_p_unw/(this_alpha/this_p_unw)<1)) {
	    opt_alpha = this_alpha;
	    opt_t_trial = this_t_trial;
	    opt_p_unw = this_p_unw;
	    opt_p_eff = this_p_eff;
	    opt_p_eff_sign = this_p_eff_sign;
	  }
	}
	sum_t_trial += opt_t_trial;
	sum_p_unw += opt_p_unw;
	sum_p_eff += opt_p_eff;
	sum_p_eff_sign += opt_p_eff_sign;
	sum_alpha += opt_alpha;
      }
      mean_manual_alpha_fscan[fi] = sum_alpha/sum_p_unw*pow(sum_p_eff_sign/sum_p_eff,2);
      mean_manual_events_fscan[fi] = 60*60*24/(sum_t_trial/sum_p_unw);
      mean_manual_eff_events_fscan[fi] = mean_manual_events_fscan[fi]*mean_manual_alpha_fscan[fi];
    }








    
    //make nice final table IV (for chosen emax)
    if (timing_statistics>3) {
      std::cout << "┌──────────────────────────────┬──────────────┬─────────────────────────────────────────────────────────┐" << std::endl;
      std::cout << "│    sampling contribution     │              │                                                         │" << std::endl;
      std::cout << "│ xsec*h  efficiency  stat.dil │ time |xsec*h|│ subprocess                                              │" << std::endl;
      std::cout << "├──────────────────────────────┼──────────────┼─────────────────────────────────────────────────────────┤" << std::endl;
      const auto default_precision{std::cout.precision()};
      std::cout << std::setprecision(3);
      for (auto const& [key, val] : time_map) {
        if (key.rfind("sum_PS_", 0) != 0) continue;
        std::string sub_name = key.substr(7);
        std::cout<<"│ "<< std::left<<std::setw(9) << xsec_map[sub_name]<<"";
        std::cout<<" "<< std::left<<std::setw(9) << chosen_efficiency_map[sub_name]<<"";
        std::cout<<" "<< std::left<<std::setw(8) << chosen_alpha_map[sub_name]<<" ";
        double this_sepsum = time_map["sum_overhead_after_"+sub_name]+time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_total_"+sub_name];
        std::cout <<"│ " <<std::right<<std::setw(2) << round(this_sepsum/(total["overhead_after"]+total["total"])*100)<<"%  ";
        std::cout <<" " <<std::right<<std::setw(4) << round(dabs(xsec_map[sub_name])/plain_xsec_sum*1000.)/10.<<"%  ";
        std::cout <<"│ " <<std::left<<std::setw(55)<< sub_name << " │" << std::endl;
      }
      std::cout << std::setprecision(default_precision);
      std::cout << "└──────────────────────────────┴──────────────┴─────────────────────────────────────────────────────────┘" << std::endl << std::endl;
    }

    //make nice final table III (for chosen emax)
    if (timing_statistics>2) {
      std::cout << "┌─────────────────────────┬──────────────┬─────────────────────────────────────────────────────────┐" << std::endl;
      std::cout << "│    time contribution    │              │                                                         │" << std::endl;
      std::cout << "│ PS   ME   ov.h.  shower │ time |xsec*h|│ subprocess                                              │" << std::endl;
      std::cout << "├─────────────────────────┼──────────────┼─────────────────────────────────────────────────────────┤" << std::endl;
      for (auto const& [key, val] : time_map) {
        if (key.rfind("sum_PS_", 0) != 0) continue;
        std::string sub_name = key.substr(7);
        if (number_map["n_total_"+sub_name]>0) {
          double t_ME = time_map["sum_ME_"+sub_name]/number_map["n_total_"+sub_name];
          double t_PS = time_map["sum_PS_"+sub_name]/number_map["n_total_"+sub_name];
          double t_ov_during = (time_map["sum_total_"+sub_name]-time_map["sum_ME_"+sub_name]-time_map["sum_PS_"+sub_name]-time_map["sum_PDF_"+sub_name])/number_map["n_total_"+sub_name];
          double t_ov_after = (time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_overhead_after_"+sub_name])/number_map["n_total_"+sub_name];
          double sum = t_ME+t_PS+t_ov_during+t_ov_after;
          msg_Info()<<"│ "<< std::right<<std::setw(2) << round(t_PS/sum*100)<<"% ";
          msg_Info()<<" "<< std::right<<std::setw(2) << round(t_ME/sum*100)<<"% ";
          msg_Info()<<"  "<< std::right<<std::setw(2) << round(t_ov_during/sum*100)<<"% ";
          msg_Info()<<"   "<< std::right<<std::setw(2) << round(t_ov_after/sum*100)<<"%   ";
        } else {
          std::cout << "│ --   --    --     --   ";
        }
        double this_sepsum = time_map["sum_overhead_after_"+sub_name]+time_map["sum_overhead_after_kept_"+sub_name]+time_map["sum_total_"+sub_name];
        std::cout <<"│ " <<std::right<<std::setw(2) << round(this_sepsum/(total["overhead_after"]+total["total"])*100)<<"%  ";
        std::cout <<" " <<std::right<<std::setw(4) << round(dabs(xsec_map[sub_name])/plain_xsec_sum*1000.)/10.<<"%  ";
        std::cout <<"│ " <<std::left<<std::setw(55)<< sub_name << " │" << std::endl;
      }
      std::cout << "└─────────────────────────┴──────────────┴─────────────────────────────────────────────────────────┘" << std::endl;
    }

    //make nice final table II
    if (timing_statistics>1) {
      std::cout << "┌─────────────┬─────────────────────────────────────────────┐" << std::endl;
      std::cout << "│ Max_Epsilon │ PS time   ME time   overhead    shower etc  │" << std::endl;
      std::cout << "├─────────────┼─────────────────────────────────────────────┤" << std::endl;
      for(int i=0; i < epsilon_values.size(); i++){
        std::cout << "│ 1e" <<std::left<<std::setw(9)<< epsilon_values[i] << " │ ";
        msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["PS"]/mean_manual_events_time[i]["sum"]*100)<<"%      ";
        msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["ME"]/mean_manual_events_time[i]["sum"]*100)<<"%      ";
        msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["ov_during"]/mean_manual_events_time[i]["sum"]*100)<<"%       ";
        msg_Info()<<" "<< std::right<<std::setw(2) << round(mean_manual_events_time[i]["ov_after"]/mean_manual_events_time[i]["sum"]*100)<<"%         │" << std::endl;
      }
      std::cout << "└─────────────┴─────────────────────────────────────────────┘" << std::endl;
    }

    //make nice final table I.5
    std::cout << "┌─────────────────────────────┬─────────────────────────────────┐" << std::endl;
    std::cout << "│                             │      Corresponding region       │" << std::endl;
    std::cout << "│ Fraction       Events/day   │ Eff. events/day     Sample size │" << std::endl;
    std::cout << "├─────────────────────────────┼─────────────────────────────────┤" << std::endl;
    const auto default_precision{std::cout.precision()};
    std::cout << std::setprecision(4);
    for(int i=0; i < fraction_values.size(); i++){
      std::cout << "│ 1e" <<std::left<<std::setw(9)<< fraction_values[i] << "    " <<std::setw(12)<< mean_manual_events_fscan[i] << " │ ";
      std::cout <<std::setw(15)<< mean_manual_eff_events_fscan[i] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha_fscan[i] << " │ " << std::endl;
    }
    std::cout << "└─────────────────────────────┴─────────────────────────────────┘" << std::endl;
    
    //make nice final table I
    std::cout << "┌─────────────────────────────┬─────────────────────────────────┬─────────────────────────────────┐" << std::endl;
    std::cout << "│                             │          Average region         │       Large weight region       │" << std::endl;
    std::cout << "│                             │           (fraction: 1)         │        (fraction: "<< timing_statistics_large_weight_fraction <<")        │" << std::endl;
    std::cout << "│ Max_Epsilon    events/day   │ eff. events/day     sample size │ eff. events/day     sample size │" << std::endl;
    std::cout << "├─────────────────────────────┼─────────────────────────────────┼─────────────────────────────────┤" << std::endl;
    std::cout << "│ 0.0 (unw.) " <<std::left<< "    " <<std::setw(12)<< mean_manual_events[epsilon_values.size()] << " │ ";
    std::cout <<std::setw(15)<< mean_manual_eff_events[epsilon_values.size()] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[epsilon_values.size()] << " │ ";
    std::cout << "-                   -          " << " │" << std::endl;
    //std::cout <<std::setw(15)<< mean_manual_fraction_eff_events[epsilon_values.size()] << "     " <<std::left<<std::setw(11) << 1./mean_manual_fraction_alpha[epsilon_values.size()] << " │" << std::endl;
    for(int i=0; i < epsilon_values.size(); i++){
      std::cout << "│ 1e" <<std::left<<std::setw(9)<< epsilon_values[i] << "    " <<std::setw(12)<< mean_manual_events[i] << " │ ";
      if (i==optimal_manual_i) {
        std::cout <<std::setw(15)<< mean_manual_eff_events[i] << "<--  " <<std::left<<std::setw(11) << 1./mean_manual_alpha[i] << " │ ";
      } else {
        std::cout <<std::setw(15)<< mean_manual_eff_events[i] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[i] << " │ ";
      }
      if (i==optimal_manual_fraction_i) {
        std::cout <<std::setw(15)<< mean_manual_fraction_eff_events[i] << "<--  " <<std::left<<std::setw(11) << 1./mean_manual_fraction_alpha[i] << " │" << std::endl;
      } else {
        std::cout <<std::setw(15)<< mean_manual_fraction_eff_events[i] << "     " <<std::left<<std::setw(11) << 1./mean_manual_fraction_alpha[i] << " │" << std::endl;
      }
    }
    std::cout << "│ weighted   " <<std::left<< "    " <<std::setw(12)<< mean_manual_events[epsilon_values.size()+1] << " │ ";
    std::cout <<std::setw(15)<< mean_manual_eff_events[epsilon_values.size()+1] << "     " <<std::left<<std::setw(11) << 1./mean_manual_alpha[epsilon_values.size()+1] << " │ ";
    std::cout << "-                   -          " << " │" << std::endl;
    std::cout << "└─────────────────────────────┴─────────────────────────────────┴─────────────────────────────────┘" << std::endl;
    //show relative sudakov uncertainty for average region optimal point
    std::cout << "Relative Sudakov uncertainty for average region optimal point: +"<< (mean_manual_events_up[optimal_manual_i]/mean_manual_events[optimal_manual_i]-1)*100. <<"% -"<< (1-mean_manual_events_down[optimal_manual_i]/mean_manual_events[optimal_manual_i])*100. <<"%" << std::endl;
    std::cout << std::setprecision(default_precision);
    if (timing_statistics_det_sim!=0) {
      std::cout << "With assumed detector simulation time of "<< timing_statistics_det_sim<<" s per kept event." << std::endl;
    } else {
      std::cout << "Note: No detector simulation time was considered for the above table(s). Use TIMING_STATISTICS_DET_SIM_IN_S to set it." << std::endl;
    }

    std::map<std::string, int> whistofill_map = rpa->gen.FillsMap();
    double average_whistofill = 0;
    double xsecsum = 0;
    for (auto const& [key, val] : time_map) {
      std::string sub_name = key.substr(7);
      average_whistofill += dabs(xsec_map[sub_name])*whistofill_map[sub_name];
      xsecsum += dabs(xsec_map[sub_name]);
    }
    average_whistofill = average_whistofill/xsecsum;
    if (average_whistofill<pow(10,5)) {
      std::cout << "Warning: There are on average only "<< average_whistofill << " events generated in the integration phase per process. For a reliable prediction of the impact of Max_Epslion choises more than 1e6 are recommended. Please increase PSI: MAX_OPT during integration." << std::endl;
    } else if (average_whistofill<pow(10,6)) {
      std::cout << "Info: There are on average only "<< average_whistofill << " events generated in the integration phase per process. For a reliable prediction of the impact of Max_Epslion choises more than 1e6 are recommended. One could increase PSI: MAX_OPT during integration." << std::endl;
    }

    std::map<std::string, double> chosen_capped_fraction_map = rpa->gen.CappedFractionMap();
    double capfraction = 0;
    for (auto const& [key, val] : time_map) {
      std::string sub_name = key.substr(7);
      capfraction += dabs(xsec_map[sub_name])*chosen_capped_fraction_map[sub_name];
    }
    capfraction = capfraction/xsecsum*100;
    if (capfraction>0.01) {
      std::cout << "Warning: OVERWEIGHT_THRESHOLD is set to: " << m_ovwth << " which results in a loss of total cross section due to max weight capping of: " << std::round(capfraction*100)/100 << "%. Please consider chosing a smaller Max_Epsilon or a larger OVERWEIGHT_THRESHOLD." << std::endl;
    }
    if (generation_mode==0) {
      //there is no pilot run for weighted events and thus the corresponding overhead can not be estimated in this case. Try to convey this message in user friendly language:
      std::cout << "Warning: TIMING_STATISTICS was run in 'weighted' mode. The estimate of 'shower etc.' is only correctly determined from 'unweighted' or '(partially) unweighted' events." << std::endl;
    }
  }
  return true;
}

void Sherpa::printEffi(std::string sub_name, std::string name, std::string num_type, std::string denom_type, std::map<std::string, int> &number_map)
{
  int n_num = number_map["n_"+num_type+"_"+sub_name];
  int n_denom = number_map["n_"+denom_type+"_"+sub_name];
  double unc = GetEffiUnc(n_num, n_denom);
  msg_Info()<< sub_name << " : "<< name << " efficiency: "<<n_num/(double)n_denom<<" +-"<< unc <<std::endl;
}

double Sherpa::GetEffiUnc(int n_num, int n_denom)
{
  //from gaussion error propagation with n_denom=n_num+rest
  return pow(n_num*(n_denom-n_num)/pow(n_denom,3),0.5);
}

void Sherpa::printTime(std::string sub_name, std::string time_type, std::map<std::string, double> &time_map, std::map<std::string, int> &number_map)
{
  double w = time_map["sum_"+time_type+"_"+sub_name];
  double w2 = time_map["sum2_"+time_type+"_"+sub_name];
  int n = number_map["n_"+time_type+"_"+sub_name];
  double unc = GetUnc(w, w2, n);
  msg_Info()<< sub_name << " : "<<"Avg time per "<< time_type << " weight: "<<w/(double)n*1000<<" ms +-"<< unc*1000. <<" ms"<<std::endl;
}


double Sherpa::GetUnc(double w, double w2,int n)
{
  return pow((w2-pow(w,2)/((double)n))/((double)n*((double)n-1.0)),0.5);
}


long int Sherpa::NumberOfEvents() const
{
  return rpa->gen.NumberOfEvents();
}

const Blob_List &Sherpa::GetBlobList() const
{
  return *p_eventhandler->GetBlobs();
}

double Sherpa::GetMEWeight(const Cluster_Amplitude &ampl,const int mode) const
{
  return p_inithandler->GetMatrixElementHandler()->
    GetWeight(ampl,ATOOLS::nlo_type::lo,mode);
}

void Sherpa::DrawLogo(const bool& shouldprintversioninfo)
{
  MyStrStream version;
  version << "SHERPA v" << SHERPA_VERSION << "." << SHERPA_SUBVERSION
          << " (" << SHERPA_NAME << ")";
  msg_Info() << Frame_Header{};

  MyStrStream logo;
  logo << om::green << "                   ." << om::reset << "_";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::green << "                  .-" << om::reset << "#.";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::green << "                 .--" << om::reset << "+@.     .";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::green << "                .----" << om::reset << "@@." << om::red << "   +" << om::reset << "#-";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::green << "               .-----" << om::reset << "+@@." << om::red << " +**" << om::reset << "@-" << "         " << version.str();
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::blue  << "       :" << om::reset << "-" << om::green << "     .-------" << om::reset << "@@@" << om::red << "+***" << om::reset << "#@-";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::blue  << "      :=" << om::reset << "#*" << om::green << "   .--------" << om::reset << "+@" << om::red << "+*****" << om::reset << "@@-" << "       Monte Carlo event generator";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::blue  << "     :===" << om::reset << "@*" << om::green << " .----------" << om::red << "+******" << om::reset << "#@@-";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::blue  << "    :====" << om::reset << "#@*" << om::green << "----------" << om::red << "+********" << om::reset << "@@@-" << "     https://sherpa-team.gitlab.io";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");
  logo << om::blue  << "   :======" << om::reset << "@@*" << om::green << "--------" << om::red << "+*********" << om::reset << "#@@@-";
  msg_Info() << Frame_Line{logo.str()}; logo.str("");

  msg_Info() << Frame_Line{"                                                                            "};
  msg_Info() << Frame_Line{"         Authors:  Enrico Bothmann, Lois Flower, Christian Gutschow,        "};
  msg_Info() << Frame_Line{"                   Stefan Hoeche, Mareen Hoppe, Joshua Isaacson,            "};
  msg_Info() << Frame_Line{"                   Max Knobbe, Frank Krauss, Peter Meinzinger,              "};
  msg_Info() << Frame_Line{"                   Davide Napoletano, Alan Price, Daniel Reichelt,          "};
  msg_Info() << Frame_Line{"                   Marek Schoenherr, Steffen Schumann, Frank Siegert        "};
  msg_Info() << Frame_Line{"  Former Authors:  Gurpreet Singh Chahal, Timo Fischer, Tanju Gleisberg,    "};
  msg_Info() << Frame_Line{"                   Hendrik Hoeth, Johannes Krause, Silvan Kuttimalai,       "};
  msg_Info() << Frame_Line{"                   Ralf Kuhn, Thomas Laubrich, Sebastian Liebschner,        "};
  msg_Info() << Frame_Line{"                   Andreas Schaelicke, Holger Schulz, Jan Winter            "};
  msg_Info() << Frame_Line{"                                                                            "};
  MyStrStream citation;
  citation << "Users are kindly asked to cite " << om::bold
           << "JHEP 12 (2024) 156" << om::reset << ".";
  msg_Info() << Frame_Line{citation.str()};
  msg_Info() << Frame_Line{"                                                                            "};
  msg_Info() << Frame_Line{"This program uses a lot of genuine and original research work by others.    "};
  msg_Info() << Frame_Line{"Users are encouraged to also cite the various original publications.        "};
  msg_Info() << Frame_Line{"                                                                            "};
  msg_Info() << Frame_Footer{};
  rpa->gen.PrintGitVersion(msg->Info(), shouldprintversioninfo);
  rpa->gen.AddCitation
    (0,"The complete Sherpa package is published under \\cite{Sherpa:2024mfk}.");
}
