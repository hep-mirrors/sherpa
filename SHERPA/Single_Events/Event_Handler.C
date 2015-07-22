#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "SHERPA/Single_Events/Signal_Processes.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif
#include <unistd.h>
#include <cassert>

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/RUsage.H"


using namespace SHERPA;
using namespace ATOOLS;

static int s_retrymax(100);

Event_Handler::Event_Handler():
  m_lastparticlecounter(0), m_lastblobcounter(0), 
  m_n(0), m_addn(0), m_sum(0.0), m_sumsqr(0.0), m_maxweight(0.0),
  m_mn(0), m_msum(0.0), m_msumsqr(0.0)
{
  p_phases  = new Phase_List;
  Data_Reader reader(" ",";","!","=");
  m_checkweight = reader.GetValue<int>("CHECK_WEIGHT", 0);
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
      msg_Out()<<"WARNING in Event_Handler::AddEventPhase"
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
  msg_Out()<<"----------------------------------------------------------\n"
	    <<"-- SHERPA generates events with the following structure --\n"
	    <<"----------------------------------------------------------\n";
  msg_Out()<<"Event generation   : ";
  switch (ToType<size_t>(rpa->gen.Variable("EVENT_GENERATION_MODE"))) {
  case 0:
    msg_Out()<<"Weighted"<<std::endl;
    break;
  case 1:
    msg_Out()<<"Unweighted"<<std::endl;
    break;
  case 2:
    msg_Out()<<"Partially unweighted"<<std::endl;
    break;
  default:
    msg_Out()<<"Unknown"<<std::endl;
    break;
  }
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg_Out()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
  msg_Out()<<"---------------------------------------------------------\n";
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

bool Event_Handler::GenerateEvent(eventtype::code mode) 
{
  DEBUG_FUNC(rpa->gen.NumberOfGeneratedEvents());
  ATOOLS::ran->SaveStatus();
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
  m_blobs.push_back(p_signal);
}

bool Event_Handler::AnalyseEvent(double & weight) {
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    if ((*pit)->Type()==eph::Analysis) {
      switch ((*pit)->Treat(&m_blobs,weight)) {
      case Return_Value::Nothing :
	break;
      case Return_Value::Success : 
        Return_Value::IncCall((*pit)->Name());
	break;
      case Return_Value::Error :
        Return_Value::IncCall((*pit)->Name());
        Return_Value::IncError((*pit)->Name());
	return false;
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

int Event_Handler::IterateEventPhases(eventtype::code & mode,double & weight) {
  Phase_Iterator pit=p_phases->begin();
  int retry = 0;
  bool hardps = true;
  do {
    if ((*pit)->Type()==eph::Analysis) {
      ++pit;
      continue;
    }

    Return_Value::code rv((*pit)->Treat(&m_blobs,weight));
    if (rv!=Return_Value::Nothing)
      msg_Tracking()<<METHOD<<"(): run '"<<(*pit)->Name()<<"' -> "
		    <<rv<<std::endl;
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
  msg_Tracking()<<METHOD<<": Event ended normally.\n";
  return 0;
}

bool Event_Handler::GenerateStandardPerturbativeEvent(eventtype::code &mode)
{
  double weight = 1.;
  bool run(true);

  InitialiseSeedBlob(ATOOLS::btp::Signal_Process,
		     ATOOLS::blob_status::needs_signal);
  do {
    weight = 1.;
    switch (IterateEventPhases(mode,weight)) {
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
  double cxs((*p_signal)["Weight"]->Get<double>());
  if (!WeightIsGood(cxs)) {
    PRINT_INFO("Invalid weight w="<<cxs<<". Rejecting event.");
    return false;
  }
  m_n      += trials+m_addn;
  m_sum    += cxs;
  m_sumsqr += sqr(cxs);
  m_addn    = 0.0;

  return AnalyseEvent(weight);
}

bool Event_Handler::GenerateMinimumBiasEvent(eventtype::code & mode) {
  double weight = 1.;
  bool run(true);

  InitialiseSeedBlob(ATOOLS::btp::Soft_Collision,
		     ATOOLS::blob_status::needs_minBias);
  do {
    weight = 1.;
    switch (IterateEventPhases(mode,weight)) {
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

  double xs((*p_signal)["Weight"]->Get<double>());
  m_n++;
  m_sum    += xs;
  m_sumsqr += sqr(xs);
  msg_Tracking()<<METHOD<<" for event with xs = "<<(xs/1.e9)<<" mbarn.\n";
  return AnalyseEvent(weight);
}


bool Event_Handler::GenerateHadronDecayEvent(eventtype::code & mode) {
  double weight = 1.;
  bool run(true);

  Data_Reader read(" ",";","!","=");
  int mother_kf(0);
  if (!read.ReadFromFile(mother_kf,"DECAYER")) {
    THROW(fatal_error,"Didn't find DECAYER=<PDG_CODE> in parameters.");
  }
  Flavour mother_flav(mother_kf);
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
    weight = 1.;
    switch (IterateEventPhases(mode,weight)) {
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

  return AnalyseEvent(weight);
}

void Event_Handler::Finish() {
  if (this==NULL) return;
  MPISync();
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
  double xs(TotalXSMPI()), err(TotalErrMPI());
  std::string res;
  MyStrStream conv;
  conv<<om::bold<<"Total XS"<<om::reset<<" is "
      <<om::blue<<om::bold<<xs<<" pb"<<om::reset<<" +- ( "
      <<om::red<<err<<" pb = "<<((int(err/xs*10000))/100.0)
      <<" %"<<om::reset<<" )";
  getline(conv,res);
  int md(msg->Modifiable()?26:-4);
  msg_Out()<<om::bold<<'+'<<std::string(res.length()-md,'-')<<"+\n";
  msg_Out()<<'|'<<std::string(res.length()-md,' ')<<"|\n";
  msg_Out()<<'|'<<om::reset<<"  "<<res<<"  "<<om::bold<<"|\n";
  msg_Out()<<'|'<<std::string(res.length()-md,' ')<<"|\n";
  msg_Out()<<'+'<<std::string(res.length()-md,'-')<<'+'<<om::reset<<std::endl;
}

void Event_Handler::MPISync()
{
  m_mn=m_n;
  m_msum=m_sum;
  m_msumsqr=m_sumsqr;
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    double values[3];
    values[0]=m_mn;
    values[1]=m_msum;
    values[2]=m_msumsqr;
    mpi->MPIComm()->Allreduce(MPI_IN_PLACE,values,3,MPI::DOUBLE,MPI::SUM);
    mpi->MPIComm()->Allreduce(MPI_IN_PLACE,&m_maxweight,1,MPI::DOUBLE,MPI::MAX);
    m_mn=values[0];
    m_msum=values[1];
    m_msumsqr=values[2];
  }
#endif
  size_t currentrss=GetCurrentRSS();
  if (m_lastrss==0) m_lastrss=currentrss;
  else if (currentrss>m_lastrss+ToType<int>
      (rpa->gen.Variable("MEMLEAK_WARNING_THRESHOLD"))) {
    msg_Error()<<METHOD<<"() {\n"<<om::bold<<"  Memory usage increased by "
	       <<(currentrss-m_lastrss)/(1<<20)<<" MB,"
	       <<" now "<<currentrss/(1<<20)<<" MB.\n"
	       <<om::red<<"  This might indicate a memory leak!\n"
	       <<"  Please monitor this process closely.\n"<<om::reset<<"}"<<std::endl;
    m_lastrss=currentrss;
  }
}

double Event_Handler::TotalXS()
{
  if (m_n==0.0) return 0.0;
  return m_sum/m_n;
}


double Event_Handler::TotalVar()
{
  if (m_n<=1) return sqr(TotalXS());
  return (m_sumsqr-m_sum*m_sum/m_n)/(m_n-1);
}


double Event_Handler::TotalErr()
{
  if (m_n<=1) return TotalXS();
  if (ATOOLS::IsEqual
      (m_sumsqr*m_n,m_sum*m_sum,1.0e-6)) return 0.0;
  return sqrt((m_sumsqr-m_sum*m_sum/m_n)/(m_n-1)/m_n);
}

double Event_Handler::TotalXSMPI()
{
  MPISync();
  if (m_mn==0.0) return 0.0;
  return m_msum/m_mn;
}


double Event_Handler::TotalVarMPI()
{
  if (m_mn<=1) return sqr(TotalXSMPI());
  return (m_msumsqr-m_msum*m_msum/m_mn)/(m_mn-1);
}


double Event_Handler::TotalErrMPI()
{
  if (m_mn<=1) return TotalXS();
  if (ATOOLS::IsEqual
      (m_msumsqr*m_mn,m_msum*m_msum,1.0e-6)) return 0.0;
  return sqrt((m_msumsqr-m_msum*m_msum/m_mn)/(m_mn-1)/m_mn);
}

bool Event_Handler::WeightIsGood(const double& weight)
{
  if (IsBad(weight)) return false;

  if (m_checkweight && fabs(weight)>m_maxweight) {
    m_maxweight=fabs(weight);
    std::string ranfilename="random.dat";
    if (ATOOLS::msg->LogFile()!="") ranfilename=ATOOLS::msg->LogFile()+"."+ranfilename;
    ATOOLS::ran->WriteOutSavedStatus(ranfilename.c_str());
    std::ofstream outstream(ranfilename.c_str(), std::fstream::app);
    outstream<<std::endl;
    outstream<<"# Wrote status for weight="<<weight<<" in event "<<rpa->gen.NumberOfGeneratedEvents()+1<<std::endl;
    outstream.close();
  }

  return true;
}
