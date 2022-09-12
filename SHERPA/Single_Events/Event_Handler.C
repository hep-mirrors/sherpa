#include "SHERPA/Single_Events/Event_Handler.H"
#include "SHERPA/Main/Filter.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/RUsage.H"
#include "SHERPA/Single_Events/Signal_Processes.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif

#include <signal.h>
#include <unistd.h>

#include <cassert>

using namespace SHERPA;
using namespace ATOOLS;

static int s_retrymax(100);

Event_Handler::Event_Handler():
  m_lastparticlecounter(0), m_lastblobcounter(0),
  m_n(0), m_mn(0), m_addn(0), m_maxweight(0.0),
  m_wgtmapsum{0.0}, m_wgtmapsumsqr{0.0},
  m_mwgtmapsum{0.0}, m_mwgtmapsumsqr{0.0},
  p_filter(NULL), p_variations(NULL)
{
  p_phases  = new Phase_List;
  Settings& s = Settings::GetMainSettings();
  m_checkweight = s["CHECK_WEIGHT"].SetDefault(0).Get<int>();
  m_decayer = s["DECAYER"].SetDefault(kf_none).Get<int>();
  m_lastrss=0;
}

Event_Handler::~Event_Handler()
{
  Reset();
  m_blobs.Clear();
  EmptyEventPhases();

  if (p_phases)   { delete p_phases;   p_phases   = NULL; }
}

void Event_Handler::AddEventPhase(Event_Phase_Handler * phase)
{
  eph::code type   = phase->Type();
  std::string name = phase->Name();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    if ((type==(*pit)->Type()) && (name==(*pit)->Name())) {
      msg_Info()<<"WARNING in Event_Handler::AddEventPhase"
		<<"("<<type<<":"<<name<<") "
		<<"already included."<<std::endl;
      return;
    }
  }
  msg_Tracking()<<"Event_Handler::AddEventPhase"
		<<"("<<type<<":"<<name<<")."<<std::endl;
  p_phases->push_back(phase);
}

void Event_Handler::EmptyEventPhases()
{
  if (p_phases) {
    while (!p_phases->empty()) {
      delete p_phases->back();
      p_phases->pop_back();
    }
  }
}

void Event_Handler::PrintGenericEventStructure()
{
  if (!msg_LevelIsInfo()) return;
  msg_Info()<<"----------------------------------------------------------\n"
	    <<"-- SHERPA generates events with the following structure --\n"
	    <<"----------------------------------------------------------\n";
  msg_Info()<<"Event generation   : ";
  switch (ToType<size_t>(rpa->gen.Variable("EVENT_GENERATION_MODE"))) {
  case 0:
    msg_Info()<<"Weighted"<<std::endl;
    break;
  case 1:
    msg_Info()<<"Unweighted"<<std::endl;
    break;
  case 2:
    msg_Info()<<"Partially unweighted"<<std::endl;
    break;
  default:
    msg_Info()<<"Unknown"<<std::endl;
    break;
  }
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg_Info()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
  if (p_variations && !p_variations->GetParametersVector()->empty()) {
    msg_Info()<<"Reweighting        : "
	     <<p_variations->GetParametersVector()->size()<<" variations.\n";
  }
  msg_Info()<<"---------------------------------------------------------\n";
}

void Event_Handler::Reset()
{
  m_sblobs.Clear();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit)
    (*pit)->CleanUp();
  m_blobs.Clear();
  if (Particle::Counter()>m_lastparticlecounter || 
      Blob::Counter()>m_lastblobcounter) {
    msg_Error()<<METHOD<<"(): "<<Particle::Counter()
               <<" particles and "<<Blob::Counter()
               <<" blobs undeleted. Continuing.\n";
    m_lastparticlecounter=Particle::Counter();
    m_lastblobcounter=Blob::Counter();
  }
  Blob::Reset();
  Particle::Reset();
  Flow::ResetCounter();
}

void Event_Handler::ResetNonPerturbativePhases()
{
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    if ((*pit)->Type()==eph::Hadronization) {
      (*pit)->CleanUp();
    }
  }
}

bool Event_Handler::GenerateEvent(eventtype::code mode) 
{
  DEBUG_FUNC(rpa->gen.NumberOfGeneratedEvents());
  ATOOLS::ran->SaveStatus();
  if (m_checkweight&4 && rpa->gen.NumberOfGeneratedEvents()==0)
    WriteRNGStatus("random","");
#ifdef USING__PYTHIA
  Lund_Interface::SaveStatus();
#endif
  if (!rpa->gen.CheckTime()) {
    msg_Error()<<ATOOLS::om::bold
                     <<"\n\nEvent_Handler::GenerateEvent("<<mode<<"): "
                     <<ATOOLS::om::reset<<ATOOLS::om::red
                     <<"Timeout. Interrupt event generation."
                     <<ATOOLS::om::reset<<std::endl;
    kill(getpid(),SIGINT);
  }
  switch (mode) {
  case eventtype::StandardPerturbative:
  case eventtype::EventReader:
    return GenerateStandardPerturbativeEvent(mode);
  case eventtype::MinimumBias:
    return GenerateMinimumBiasEvent(mode);
  case eventtype::HadronDecay:
    return GenerateHadronDecayEvent(mode);
  }
  return false;
} 

void Event_Handler::InitialiseSeedBlob(ATOOLS::btp::code type,
				       ATOOLS::blob_status::code status) {
  p_signal=new Blob();
  p_signal->SetType(type);
  p_signal->SetId();
  p_signal->SetStatus(status);
  p_signal->AddData("Trials",new Blob_Data<double>(0));
  p_signal->AddData("WeightsMap",new Blob_Data<Weights_Map>({}));
  p_signal->AddData("Weight_Norm",new Blob_Data<double>(1.));
  m_blobs.push_back(p_signal);
}

bool Event_Handler::AnalyseEvent() {
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    if ((*pit)->Type()==eph::Analysis) {
      switch ((*pit)->Treat(&m_blobs)) {
      case Return_Value::Nothing :
	break;
      case Return_Value::Success : 
        Return_Value::IncCall((*pit)->Name());
	break;
      case Return_Value::Error :
        Return_Value::IncCall((*pit)->Name());
        Return_Value::IncError((*pit)->Name());
        return false;
      case Return_Value::New_Event :
	{
	  double trials=(*p_signal)["Trials"]->Get<double>();
	  const auto& lastwgtmap = (*p_signal)["WeightsMap"]->Get<Weights_Map>();
	  m_n -= trials;
	  m_addn = trials;
	  m_wgtmapsum -= lastwgtmap;
	  m_wgtmapsumsqr -= lastwgtmap*lastwgtmap;
	  Return_Value::IncCall((*pit)->Name());
	  Return_Value::IncNewEvent((*pit)->Name());
	  Reset();
	  return false;
	}
      default:
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"  Unknown return value for 'Treat',\n"
		   <<"  Will continue and hope for the best.\n";
	return false;
      }
    }
  }
  return true;
}

int Event_Handler::IterateEventPhases(eventtype::code & mode) {
  Phase_Iterator pit=p_phases->begin();
  int retry = 0;
  bool hardps = true, filter = p_filter!=NULL;
  do {
    if ((*pit)->Type()==eph::Analysis || (*pit)->Type()==eph::Userhook) {
      ++pit;
      continue;
    }
    if ((*pit)->Type()==eph::Hadronization && filter) {
      msg_Debugging()<<"Filter kicks in now: "<<m_blobs.back()->Type()<<".\n";
      if ((*p_filter)(&m_blobs)) {
	msg_Debugging()<<METHOD<<": filters accepts event.\n";
	filter = false;
      }
      else {
	msg_Debugging()<<METHOD<<": filter rejects event.\n";
	Return_Value::IncNewEvent("Filter");
	if (p_signal) m_addn+=(*p_signal)["Trials"]->Get<double>();
	Reset();
	return 2;
      }
    }
    DEBUG_INFO("Treating "<<(*pit)->Name());
    Return_Value::code rv((*pit)->Treat(&m_blobs));
    if (rv!=Return_Value::Nothing)
      msg_Tracking()<<METHOD<<"(): run '"<<(*pit)->Name()<<"' -> "
                    <<rv<<std::endl;
      msg_Debugging()<<" -> "<<rv<<" ("<<m_blobs.size()<<" blobs)"<<std::endl;
    switch (rv) {
    case Return_Value::Success : 
      if (mode==eventtype::StandardPerturbative &&
	  (*pit)->Name().find("Jet_Evolution")==0 && hardps) {
	m_sblobs.Clear();
	m_sblobs=m_blobs.Copy();
	hardps=false;
      }
      Return_Value::IncCall((*pit)->Name());
      msg_Debugging()<<m_blobs;
      pit=p_phases->begin();
      break;
    case Return_Value::Nothing :
      ++pit;
      break;
    case Return_Value::Retry_Phase :
      Return_Value::IncCall((*pit)->Name());
      Return_Value::IncRetryPhase((*pit)->Name());
      break;
    case Return_Value::Retry_Event :
      if (retry <= s_retrymax) {
        retry++;
        Return_Value::IncCall((*pit)->Name());
        Return_Value::IncRetryEvent((*pit)->Name());
        ResetNonPerturbativePhases();
	if (mode==eventtype::StandardPerturbative) {
	  m_blobs.Clear();
	  m_blobs=m_sblobs.Copy();
	  p_signal=m_blobs.FindFirst(btp::Signal_Process);
	  if (p_signal) {
	    pit=p_phases->begin();
	    break;
	  }
	}
      }
      else {
	msg_Error()<<METHOD<<"(): No success after "<<s_retrymax
		   <<" trials. Request new event.\n";
      }
    case Return_Value::New_Event :
      Return_Value::IncCall((*pit)->Name());
      Return_Value::IncNewEvent((*pit)->Name());
      if (p_signal) m_addn+=(*p_signal)["Trials"]->Get<double>();
      Reset();
      return 2;
    case Return_Value::Error :
      Return_Value::IncCall((*pit)->Name());
      Return_Value::IncError((*pit)->Name());
      return 3;
    default:
      THROW(fatal_error,"Invalid return value");
    }
  } while (pit!=p_phases->end());
  msg_Tracking()<<METHOD<<": Event phases ended normally.\n";

  msg_Tracking()<<METHOD<<": Running user hooks now.\n";
  for (size_t i=0; i<p_phases->size(); ++i) {
    Event_Phase_Handler* phase=(*p_phases)[i];
    if (phase->Type()!=eph::Userhook) continue;

    Return_Value::code rv(phase->Treat(&m_blobs));
    if (rv!=Return_Value::Nothing)
      msg_Tracking()<<METHOD<<"(): ran '"<<phase->Name()<<"' -> "
		    <<rv<<std::endl;
    switch (rv) {
    case Return_Value::Success :
      Return_Value::IncCall(phase->Name());
      msg_Debugging()<<m_blobs;
      break;
    case Return_Value::Nothing :
      break;
    case Return_Value::New_Event :
      Return_Value::IncCall(phase->Name());
      Return_Value::IncNewEvent(phase->Name());
      if (p_signal) m_addn+=(*p_signal)["Trials"]->Get<double>();
      Reset();
      return 2;
    case Return_Value::Error :
      Return_Value::IncCall(phase->Name());
      Return_Value::IncError(phase->Name());
      return 3;
    default:
      THROW(fatal_error,"Invalid return value");
    }
  }
  msg_Tracking()<<METHOD<<": User hooks ended normally.\n";

  return 0;
}

bool Event_Handler::GenerateStandardPerturbativeEvent(eventtype::code &mode)
{
  DEBUG_FUNC(mode);
  bool run(true);

  InitialiseSeedBlob(ATOOLS::btp::Signal_Process,
		     ATOOLS::blob_status::needs_signal);
  do {
    switch (IterateEventPhases(mode)) {
    case 3:
      return false;
    case 2:
      InitialiseSeedBlob(ATOOLS::btp::Signal_Process,
			 ATOOLS::blob_status::needs_signal);
      break;
    case 1:
      m_blobs.Clear(p_signal);
      p_signal->SetStatus(blob_status::internal_flag | 
			  blob_status::needs_signal);
      break;
    case 0:
      run = false;
      break;
    }
  } while (run);

  if (mode==eventtype::EventReader) {
    if (p_signal->NOutP()==0) return false;
  }
  else {
    if (!m_blobs.FourMomentumConservation()) {
      msg_Debugging()<<m_blobs<<"\n";
      msg_Error()<<METHOD<<"(): "
		 <<"Four momentum not conserved. Rejecting event.\n";
      return false;
    }
  }

  double trials((*p_signal)["Trials"]->Get<double>());
  p_signal->AddData("Trials",new Blob_Data<double>(trials+m_addn));

  Weights_Map wgtmap((*p_signal)["WeightsMap"]->Get<Weights_Map>());

  if (!WeightsAreGood(wgtmap)) {
    PRINT_INFO("Invalid weight w="<<wgtmap.Nominal()<<". Rejecting event.");
    return false;
  }
  m_n      += trials+m_addn;
  m_addn    = 0.0;
  m_wgtmapsum += wgtmap;
  m_wgtmapsumsqr += wgtmap*wgtmap;

  return AnalyseEvent();
}

bool Event_Handler::GenerateMinimumBiasEvent(eventtype::code & mode) {
  bool run(true);
  InitialiseSeedBlob(ATOOLS::btp::Soft_Collision,
		     ATOOLS::blob_status::needs_minBias);
  do {
    switch (IterateEventPhases(mode)) {
    case 3:
      return false;
    case 2:
    case 1:
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
        (*pit)->CleanUp();
      }
      m_blobs.Clear();
      if (Particle::Counter()>m_lastparticlecounter || 
	  Blob::Counter()>m_lastblobcounter) {
	msg_Error()<<METHOD<<"(): "<<Particle::Counter()
		   <<" particles and "<<Blob::Counter()
		   <<" blobs undeleted. Continuing.\n";
	m_lastparticlecounter=Particle::Counter();
	m_lastblobcounter=Blob::Counter();
      }
      InitialiseSeedBlob(ATOOLS::btp::Soft_Collision,
			 ATOOLS::blob_status::needs_minBias);
      break;
    case 0:
      run = false;
      break;
    }
  } while (run);

  Weights_Map wgtmap((*p_signal)["WeightsMap"]->Get<Weights_Map>());
  m_n++;
  m_wgtmapsum += wgtmap;
  m_wgtmapsumsqr += wgtmap*wgtmap;

  return AnalyseEvent();
}


bool Event_Handler::GenerateHadronDecayEvent(eventtype::code & mode) {
  bool run(true);
  if (m_decayer == kf_none) {
    THROW(fatal_error,"Didn't find DECAYER=<PDG_CODE> in parameters.");
  }
  Flavour mother_flav(m_decayer);
  mother_flav.SetStable(false);
  rpa->gen.SetEcms(mother_flav.HadMass());

  InitialiseSeedBlob(ATOOLS::btp::Hadron_Decay,
                     ATOOLS::blob_status::needs_hadrondecays);
  Vec4D mom(mother_flav.HadMass(), 0., 0., 0.);
  Particle* mother_in_part=new Particle(1, mother_flav, mom);
  Particle* mother_part=new Particle(1, mother_flav, mom);
  mother_part->SetTime();
  mother_part->SetFinalMass(mother_flav.HadMass());
  mother_in_part->SetStatus(part_status::decayed);
  p_signal->SetStatus(blob_status::needs_hadrondecays);
  p_signal->AddToInParticles(mother_in_part);
  p_signal->AddToOutParticles(mother_part);
  
  do {
    switch (IterateEventPhases(mode)) {
    case 3:
      return false;
    case 2:
      InitialiseSeedBlob(ATOOLS::btp::Hadron_Decay,
                         ATOOLS::blob_status::needs_hadrondecays);
      mother_in_part=new Particle(1, mother_flav, mom);
      mother_part=new Particle(1, mother_flav, mom);
      mother_part->SetTime();
      mother_part->SetFinalMass(mother_flav.HadMass());
      mother_in_part->SetStatus(part_status::decayed);
      p_signal->SetStatus(blob_status::needs_hadrondecays);
      p_signal->AddToInParticles(mother_in_part);
      p_signal->AddToOutParticles(mother_part);
      break;
    case 1:
      m_blobs.Clear(p_signal);
      p_signal->SetStatus(blob_status::internal_flag | 
                          blob_status::needs_minBias);
      break;
    case 0:
      run = false;
      break;
    }
  } while (run);

  return AnalyseEvent();
}

void Event_Handler::Finish() {
  msg_Info()<<"In Event_Handler::Finish : "
	    <<"Summarizing the run may take some time.\n";
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    (*pit)->Finish(std::string("Results"));
    (*pit)->CleanUp();
  }
  m_blobs.Clear();
  m_sblobs.Clear();
  if (Particle::Counter()>m_lastparticlecounter || 
      Blob::Counter()>m_lastblobcounter) {
    msg_Error()<<"ERROR in "<<METHOD<<":\n"
	       <<"   After event : "<<Particle::Counter()
	       <<" / "<<Blob::Counter()
	       <<" particles / blobs undeleted !\n"
	       <<"   Continue and hope for the best.\n";
    m_lastparticlecounter=Particle::Counter();
    m_lastblobcounter=Blob::Counter();
  }
  Blob::Reset();
  // Obtain absolute (variation) weights.
  MPISyncXSAndErrMaps();
  Weights_Map xs_wgtmap = TotalXSMPI();
  Weights_Map err_wgtmap = TotalErrMPI();
  std::map<std::string, double> xs_wgts;
  xs_wgtmap.FillVariations(xs_wgts);
  std::map<std::string, double> err_wgts;
  err_wgtmap.FillVariations(err_wgts);
  if (Settings::GetMainSettings()["OUTPUT_ME_ONLY_VARIATIONS"]
	  .SetDefault(true)
	  .Get<bool>()) {
    xs_wgtmap.FillVariations(xs_wgts, Variations_Source::main);
    err_wgtmap.FillVariations(err_wgts, Variations_Source::main);
  }

  // Find longest weights name
  static const std::string name_column_title {"Nominal or variation name"};
  size_t max_weight_name_size {name_column_title.size()};
  for (const auto& kv : xs_wgts)
    max_weight_name_size = std::max(max_weight_name_size, kv.first.size());

  // Calculate columns widths
  const size_t xs_size {12};
  const size_t reldev_size {12};
  const size_t abserr_size {13};
  const size_t relerr_size {12};
  const size_t table_size {max_weight_name_size + xs_size + reldev_size +
			   abserr_size + relerr_size};

  // Print cross section table header.
  msg_Out() << std::string(table_size, '-') << '\n';
  msg_Out() << std::left << std::setw(max_weight_name_size)
	    << "Nominal or variation name";
  msg_Out() << std::right << std::setw(12) << "XS [pb]";
  msg_Out() << std::right << std::setw(12) << "RelDev";
  msg_Out() << std::right << std::setw(13) << "AbsErr [pb]";
  msg_Out() << std::right << std::setw(12) << "RelErr" << '\n';
  msg_Out() << std::string(table_size, '-') << '\n';
  // Define table row printer.
  auto printxs = [max_weight_name_size, xs_size, reldev_size, abserr_size,
		  relerr_size](const std::string& name, double xs, double nom,
			       double err) {
    msg_Out() << om::bold << std::left << std::setw(max_weight_name_size)
	      << name << om::reset << std::right << om::blue << om::bold
	      << std::setw(xs_size) << xs << om::reset << om::brown
	      << std::setw(reldev_size - 2)
	      << ((int((xs - nom) / nom * 10000)) / 100.0) << " %" << om::red
	      << std::setw(abserr_size) << err << std::setw(relerr_size - 2)
	      << ((int(err / xs * 10000)) / 100.0) << " %" << om::reset
	      << std::endl;
  };

  // Print nominal cross section and variations.
  double nom = xs_wgtmap.Nominal();
  printxs("Nominal", nom, nom, err_wgtmap.Nominal());
  for (const auto& kv : xs_wgts) {
    const double xs = kv.second;
    const double err = err_wgts[kv.first];
    printxs(kv.first, xs, nom, err);
  }

  // Print cross section table footer.
  msg_Out()<<std::string(table_size,'-')<<'\n';
}

void Event_Handler::MPISync()
{
  m_mn = m_n;
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(&m_mn, 1, MPI_DOUBLE, MPI_SUM);
    if (!(m_checkweight & 2))
      mpi->Allreduce(&m_maxweight, 1, MPI_DOUBLE, MPI_MAX);
  }
#endif
}

void Event_Handler::MPISyncXSAndErrMaps()
{
  MPISync();
  m_mwgtmapsum = m_wgtmapsum;
  m_mwgtmapsumsqr = m_wgtmapsumsqr;
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    m_mwgtmapsum.MPI_Allreduce();
    m_mwgtmapsumsqr.MPI_Allreduce();
  }
#endif
}

void Event_Handler::PerformMemoryMonitoring()
{
  size_t currentrss=GetCurrentRSS();
  if (m_lastrss==0) m_lastrss=currentrss;
  else if (currentrss>m_lastrss+ToType<int>
      (rpa->gen.Variable("MEMLEAK_WARNING_THRESHOLD"))) {
    msg_Error()<<"\n"<<om::bold<<"    Memory usage increased by "
	       <<(currentrss-m_lastrss)/(1<<20)<<" MB,"
	       <<" now "<<currentrss/(1<<20)<<" MB.\n"
	       <<om::red<<"    This might indicate a memory leak!\n"
	       <<"    Please monitor this process closely."<<om::reset<<std::endl;
    m_lastrss=currentrss;
  }
}

Uncertain<double> Event_Handler::TotalNominalXS()
{
  if (m_n == 0.0)
    return {0.0, 0.0};

  const double sum_nominal {m_wgtmapsum.Nominal()};
  const double xs {sum_nominal / m_n};
  if (m_n <= 1)
    return {xs, xs};

  const double sumsqr_nominal {m_wgtmapsumsqr.Nominal()};
  if (ATOOLS::IsEqual(sumsqr_nominal * m_n, sum_nominal * sum_nominal, 1.0e-6))
    return {xs, 0.0};

  const double numerator {sumsqr_nominal * m_n - sum_nominal * sum_nominal};
  return {xs, sqrt(numerator / (m_n - 1) / (m_n * m_n))};
}

Uncertain<double> Event_Handler::TotalNominalXSMPI()
{
  MPISync();
  if (m_mn == 0.0)
    return {0.0, 0.0};

  double sum_nominal {m_wgtmapsum.Nominal()};
#ifdef USING__MPI
  if (mpi->Size() > 1)
    mpi->Allreduce(&sum_nominal, 1, MPI_DOUBLE, MPI_SUM);
#endif
  const double xs {sum_nominal / m_mn};
  if (m_mn <= 1)
    return {xs, xs};

  double sumsqr_nominal {m_wgtmapsumsqr.Nominal()};
#ifdef USING__MPI
  if (mpi->Size() > 1)
    mpi->Allreduce(&sumsqr_nominal, 1, MPI_DOUBLE, MPI_SUM);
#endif
  if (ATOOLS::IsEqual(sumsqr_nominal * m_mn, sum_nominal * sum_nominal, 1.0e-6))
    return {xs, 0.0};

  const double numerator {sumsqr_nominal * m_mn - sum_nominal * sum_nominal};
  return {xs, sqrt(numerator / (m_mn - 1) / (m_mn * m_mn))};
}

Weights_Map Event_Handler::TotalXS()
{
  if (m_n==0.0) return 0.0;
  return m_wgtmapsum/m_n;
}

Weights_Map Event_Handler::TotalErr()
{
  if (m_n<=1) return TotalXS();
  auto numerator = m_wgtmapsumsqr*m_n - m_wgtmapsum*m_wgtmapsum;
  numerator.SetZeroIfCloseToZero(1.0e-6);
  return sqrt(numerator/(m_n-1)/(m_n*m_n));
}

Weights_Map Event_Handler::TotalXSMPI()
{
  if (m_mn==0.0) return 0.0;
  return m_mwgtmapsum/m_mn;
}

Weights_Map Event_Handler::TotalErrMPI()
{
  if (m_mn<=1) return TotalXSMPI();
  auto numerator = m_mwgtmapsumsqr*m_mn - m_mwgtmapsum*m_mwgtmapsum;
  numerator.SetZeroIfCloseToZero(1.0e-6);
  return sqrt(numerator/(m_mn-1)/(m_mn*m_mn));
}

void Event_Handler::WriteRNGStatus
(const std::string &file,const std::string &message) const
{
  std::string ranfilename=file+".dat";
  if (m_checkweight&2) ranfilename=file+"."+rpa->gen.Variable("RNG_SEED")+".dat";
  if (ATOOLS::msg->LogFile()!="") ranfilename=ATOOLS::msg->LogFile()+"."+ranfilename;
  ATOOLS::ran->WriteOutSavedStatus(ranfilename.c_str());
  std::ofstream outstream(ranfilename.c_str(), std::fstream::app);
  outstream<<"\n"<<message<<"\n";
  outstream.close();
}

bool Event_Handler::WeightsAreGood(const Weights_Map& wgtmap)
{
  const auto weight = wgtmap.Nominal();
  if (IsBad(weight)) return false;

  if (m_checkweight && fabs(weight)>m_maxweight) {
    m_maxweight=fabs(weight);
    WriteRNGStatus("maxweight","# Wrote status for weight="+ToString(weight)+
		   " in event "+ToString(rpa->gen.NumberOfGeneratedEvents()+1)+
		   " trial "+ToString(rpa->gen.NumberOfTrials()-1));
  }
  if (m_checkweight & 8) {
    for (auto type : s_variations->ManagedVariationTypes()) {
      auto weights = wgtmap.Combine(type);
      const auto num_variations = s_variations->Size(type);
      for (auto i = 0; i < num_variations; ++i) {
        const auto varweight = weights.Variation(i);
        const std::string& name = s_variations->Parameters(i).Name();
        if (m_maxweights.find(name) == m_maxweights.end()) {
          m_maxweights[name] = 0.0;
        }
        if (fabs(varweight) > m_maxweights[name]) {
          m_maxweights[name] = fabs(varweight);
          WriteRNGStatus("maxweight." + name,
              "# Wrote status for weight=" + ToString(varweight) +
              " in event " +
              ToString(rpa->gen.NumberOfGeneratedEvents() + 1) +
              " trial " + ToString(rpa->gen.NumberOfTrials() - 1));
        }
      }
    }
  }

  return true;
}

std::string Event_Handler::CurrentProcess() const
{
  if (p_signal) return p_signal->TypeSpec();
  return "<unknown>";
}

