#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "SHERPA/Single_Events/Signal_Processes.H"
#ifdef USING__PYTHIA
#include "SHERPA/LundTools/Lund_Interface.H"
#endif
#include <unistd.h>
#include <cassert>
#ifdef USING__MPI
#include "mpi.h"
#endif

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Data_Reader.H"


using namespace SHERPA;
using namespace ATOOLS;

static int s_retrymax(100);

Event_Handler::Event_Handler():
  m_lastparticlecounter(0), m_lastblobcounter(0), 
  m_n(0), m_addn(0), m_sum(0.0), m_sumsqr(0.0)
#ifdef USING__MPI
  , m_mn(0), m_msum(0.0), m_msumsqr(0.0)
#endif
{
  p_phases  = new Phase_List;
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
  msg_Out()<<"----------------------------------------------------------"<<std::endl
	    <<"-- SHERPA generates events with the following structure --"<<std::endl
	    <<"----------------------------------------------------------"<<std::endl;
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg_Out()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
  msg_Out()<<"---------------------------------------------------------"<<std::endl;
}

void Event_Handler::Reset(const size_t & mode)
{
  m_sblobs.Clear();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit)
    (*pit)->CleanUp();
  m_blobs.Clear();
  if (mode==1) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Name()==std::string("Beam_Remnants")) {
	(*pit)->CleanUp(1);    
      }
    }
  }
  if (Particle::Counter()>m_lastparticlecounter || 
      Blob::Counter()>m_lastblobcounter) {
    msg_Error()<<METHOD<<"(): "<<Particle::Counter()
               <<" particles and "<<Blob::Counter()
               <<" blobs undeleted. Continuing."<<std::endl;
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
    //msg_Out()<<METHOD<<" for StandardPerturbative.\n";
    return GenerateStandardPerturbativeEvent(mode);
  case eventtype::MinimumBias:
    //msg_Out()<<METHOD<<" for MinimumBias.\n";
    return GenerateMinimumBiasEvent();
  case eventtype::HadronDecay:
    //msg_Out()<<METHOD<<" for HadronDecay.\n";
    return GenerateHadronDecayEvent();
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
	msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		   <<"  Unknown return value for 'Treat',"<<std::endl
		   <<"  Will continue and hope for the best."<<std::endl;
	return false;
      }
    }
  }
  return true;
}

int Event_Handler::IterateEventPhases(double & weight) {
  Phase_Iterator pit=p_phases->begin();
  int retry = 0;
  do {
/*    msg_Out()<<"\n";
    PRINT_VAR((*pit)->Type());
    PRINT_VAR((*pit)->Name());
    msg_Out()<<m_blobs<<"\n";*/
    if ((*pit)->Type()==eph::Analysis) {
      ++pit;
      continue;
    }

    Return_Value::code rv((*pit)->Treat(&m_blobs,weight));
//     PRINT_VAR(rv);
    if (rv!=Return_Value::Nothing)
      msg_Tracking()<<METHOD<<"(): run '"<<(*pit)->Name()<<"' -> "<<rv<<std::endl;
    switch (rv) {
    case Return_Value::Success : 
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
        Blob::Reset();
        Particle::Reset();
        Flow::ResetCounter();
        return 1;
      }
      else {
        msg_Error()<<"Too many retrials for event, generating new one.\n";
      }
    case Return_Value::New_Event : 
      Return_Value::IncCall((*pit)->Name());
      Return_Value::IncNewEvent((*pit)->Name());
      if (p_signal) m_addn+=(*p_signal)["Trials"]->Get<double>();
//       msg_Out()<<"Got New_Event, is this causing the problem?"<<std::endl;
      Reset(1);
      return 2;
    case Return_Value::Error :
      Return_Value::IncCall((*pit)->Name());
      Return_Value::IncError((*pit)->Name());
      return 3;
    default:
      THROW(fatal_error,"Invalid return value");
    }
  } while (pit!=p_phases->end());
  msg_Tracking()<<METHOD<<": Event ended normally."<<std::endl;
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
    switch (IterateEventPhases(weight)) {
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
      msg_Error()<<METHOD<<"(): Four momentum not conserved. Rejecting event."<<std::endl;
      return false;
    }
  }

  double trials((*p_signal)["Trials"]->Get<double>());
  p_signal->AddData("Trials",new Blob_Data<double>(trials+m_addn));
  double cxs((*p_signal)["Weight"]->Get<double>());
  m_n      += trials+m_addn;
  m_sum    += cxs;
  m_sumsqr += sqr(cxs);
#ifdef USING__MPI
  m_mn      += trials+m_addn;
  m_msum    += cxs;
  m_msumsqr += sqr(cxs);
#endif
  m_addn    = 0.0;

  return AnalyseEvent(weight);
}

bool Event_Handler::GenerateMinimumBiasEvent() {
  double weight = 1.;
  bool run(true);

  InitialiseSeedBlob(ATOOLS::btp::Soft_Collision,
		     ATOOLS::blob_status::needs_minBias);
  //msg_Out()<<METHOD<<" for: "<<(*p_signal)<<"\n";
  do {
    weight = 1.;
    switch (IterateEventPhases(weight)) {
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
		   <<" blobs undeleted. Continuing."<<std::endl;
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

  //msg_Out()<<METHOD<<" done with event:\n"<<m_blobs<<"\n";
  double xs((*p_signal)["Weight"]->Get<double>());
  m_n++;
  m_sum    += xs;
  m_sumsqr += sqr(xs);
#ifdef USING__MPI
  m_mn++;
  m_msum    += xs;
  m_msumsqr += sqr(xs);
#endif
  msg_Tracking()<<METHOD<<" for event with xs = "<<(xs/1.e9)<<" mbarn.\n";
  return AnalyseEvent(weight);
}


bool Event_Handler::GenerateHadronDecayEvent() {
  double weight = 1.;
  bool run(true);

  Data_Reader read(" ",";","!","=");
  int mother_kf(0);
  if (!read.ReadFromFile(mother_kf,"DECAYER")) {
    THROW(fatal_error,"Didn't find DECAYER=<PDG_CODE> in parameters.");
  }
  Flavour mother_flav;
  mother_flav.FromHepEvt(mother_kf);
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
    switch (IterateEventPhases(weight)) {
    case 3:
      return false;
    case 2:
      InitialiseSeedBlob(ATOOLS::btp::Soft_Collision,
                         ATOOLS::blob_status::needs_minBias);
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
  msg_Info()<<"In Event_Handler::Finish : Summarizing the run may take some time."<<std::endl;
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    (*pit)->Finish(std::string("Results"));
    (*pit)->CleanUp();
  }
  m_blobs.Clear();
  m_sblobs.Clear();
  if (Particle::Counter()>m_lastparticlecounter || 
      Blob::Counter()>m_lastblobcounter) {
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   After event : "<<Particle::Counter()<<" / "<<Blob::Counter()
	       <<" particles / blobs undeleted !"<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    m_lastparticlecounter=Particle::Counter();
    m_lastblobcounter=Blob::Counter();
  }
  Blob::Reset();
  double xs(TotalStatsMPI()[0]), err(TotalStatsMPI()[2]);
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
#ifdef USING__MPI
  int size=MPI::COMM_WORLD.Get_size();
  if (size>1) {
    int rank=MPI::COMM_WORLD.Get_rank();
    int cn=3;
    double *values = new double[cn];
    if (rank==0) {
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Recv(values,cn,MPI::DOUBLE,MPI::ANY_SOURCE,tag);
        m_mn+=values[0];
        m_msum+=values[1];
        m_msumsqr+=values[2];
      }
      values[0]=m_mn;
      values[1]=m_msum;
      values[2]=m_msumsqr;
      for (int tag=1;tag<size;++tag) {
	if (!exh->MPIStat(tag)) continue;
	MPI::COMM_WORLD.Send(values,cn,MPI::DOUBLE,tag,size+tag);
      }
    }
    else {
      values[0]=m_mn;
      values[1]=m_msum;
      values[2]=m_msumsqr;
      MPI::COMM_WORLD.Send(values,cn,MPI::DOUBLE,0,rank);
      MPI::COMM_WORLD.Recv(values,cn,MPI::DOUBLE,0,size+rank);
      m_mn=0.0;
      m_msum=0.0;
      m_msumsqr=0.0;
    }
    delete [] values;
  }
#endif
}

double Event_Handler::TotalXS() const
{
  return m_sum/m_n;
}


double Event_Handler::TotalVar() const
{
  return (m_sumsqr-m_sum*m_sum/m_n)/(m_n-1);
}


double Event_Handler::TotalErr() const
{
  if (m_n<=1) return TotalXS();
  else if (ATOOLS::IsEqual
           (m_sumsqr*m_n,m_sum*m_sum,1.0e-6)) return 0.0;
  else return sqrt(TotalVar()/m_n);
}

std::vector<double> Event_Handler::TotalStatsMPI()
{
  std::vector<double> stats(3, 0.0);
#ifdef USING__MPI
  MPISync();
  stats[0] = m_msum/m_mn;
  stats[1] = (m_msumsqr-m_msum*m_msum/m_mn)/(m_mn-1);
  if (m_mn<=1) stats[2] = stats[0];
  else if (IsEqual(m_msumsqr*m_mn,m_msum*m_msum,1.0e-6)) stats[2] = 0.0;
  else stats[2] = sqrt(stats[1]/m_mn);
#else
  stats[0] = TotalXS();
  stats[1] = TotalVar();
  stats[2] = TotalErr();
#endif
  return stats;
}
