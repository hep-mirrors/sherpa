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

#include "ATOOLS/Math/Random.H"

#ifdef PROFILE__Event_Handler
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;

static int s_retrymax(100);

Event_Handler::Event_Handler():
  m_lastparticlecounter(0), m_lastblobcounter(0), 
  m_n(0), m_addn(0), m_sum(0.0), m_sumsqr(0.0)
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
      msg_Out()<<"WARNING in Event_Handler::AddEventPhase("<<type<<":"<<name<<") "
		  <<"already included."<<std::endl;
      return;
    }
  }
  msg_Tracking()<<"Event_Handler::AddEventPhase("<<type<<":"<<name<<")."<<std::endl;
  p_phases->push_back(phase);
}

Event_Phase_Handler * Event_Handler::GetEventPhase(const size_t i) {
  if (i<p_phases->size()) {
    size_t count=i;
    for (Phase_Iterator pit=p_phases->begin();pit<p_phases->end();pit++) {
      if (count==0) return (*pit);
      count--;
    }
  }
  msg_Error()<<"Error in Event_Handler::GetEventPhase("<<i<<")"<<std::endl
	     <<"   Out of bounds, only "<<p_phases->size()<<" event phases."<<std::endl
	     <<"   Will return NULL."<<std::endl;
  return NULL;
}

size_t Event_Handler::NumberOfEventPhases() { return p_phases->size(); }

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

void Event_Handler::Reset()
{
  m_sblobs.Clear();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) 
    if ((*pit)->Type()!=eph::Perturbative ||
	(*pit)->Name().find("Signal_Processes")==std::string::npos) 
      (*pit)->CleanUp();
  m_blobs.Clear();
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

bool Event_Handler::GenerateEvent(int mode) 
{
  DEBUG_FUNC(rpa.gen.NumberOfDicedEvents());
  ATOOLS::ran.SaveStatus();
#ifdef USING__PYTHIA
  Lund_Interface::SaveStatus();
#endif
  if (!rpa.gen.CheckTime()) {
    msg_Error()<<ATOOLS::om::bold
                     <<"\n\nEvent_Handler::GenerateEvent("<<mode<<"): "
                     <<ATOOLS::om::reset<<ATOOLS::om::red
                     <<"Timeout. Interrupt event generation."
                     <<ATOOLS::om::reset<<std::endl;
    kill(getpid(),SIGINT);
  }
  double weight = 1.;
  switch (mode) {
  case 0:
    if (m_blobs.size()==1) {
      p_signal=m_blobs.back();
    }
    else {
      p_signal=new Blob();
      p_signal->SetType(btp::Signal_Process);
      p_signal->SetId();
      p_signal->SetStatus(blob_status::needs_signal);
      m_blobs.push_back(p_signal);
    }
    do {
      int retry(0);
      bool hardps(true);
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();) {
	if ((*pit)->Type()==eph::Analysis) {
	  ++pit;
	  continue;
	}
	msg_Debugging()<<METHOD<<"(): run '"<<(*pit)->Name()<<"'"<<std::endl;
	switch ((*pit)->Treat(&m_blobs,weight)) {
	case Return_Value::Nothing :
          DEBUG_INFO("Nothing");
	  msg_Debugging()<<m_blobs;
	  ++pit;
	  break;
	case Return_Value::Success : 
	  if ((*pit)->Name().find("Jet_Evolution")==0 && hardps) {
	    m_sblobs.Clear();
	    m_sblobs=m_blobs.Copy();
	    hardps=false;
	  }
	  rvalue.IncCall((*pit)->Name());
	  pit=p_phases->begin();
          DEBUG_INFO("Success");
          DEBUG_VAR(m_blobs.FourMomentumConservation());
	  msg_Debugging()<<m_blobs;
	  break;
	case Return_Value::Error :
          DEBUG_INFO("Error");
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncError((*pit)->Name());
	  return false;
	case Return_Value::Retry_Phase : 
          DEBUG_INFO("Retry_Phase");
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncRetryPhase((*pit)->Name());
	  pit=p_phases->begin();
	  break;
	case Return_Value::Retry_Event : 
          DEBUG_INFO("Retry_Event");
	  if (++retry<s_retrymax) {
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncRetryEvent((*pit)->Name());
          m_blobs.Clear();
	  m_blobs=m_sblobs.Copy();
	  p_signal=m_blobs.FindFirst(btp::Signal_Process);
	  pit=p_phases->begin();
	  break;
	  }
	  else {
	    msg_Error()<<METHOD<<"(): No success after "<<s_retrymax
		       <<" trials. Request new event."<<std::endl;
	  }
	case Return_Value::New_Event :
          DEBUG_INFO("New_Event");
	  {
	    Blob *sp(m_blobs.FindFirst(btp::Signal_Process));
            if (sp && (*sp)["Trials"]) m_addn+=(*sp)["Trials"]->Get<double>();
	  }
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncNewEvent((*pit)->Name());
	  Reset();
          p_signal=new Blob();
          p_signal->SetType(btp::Signal_Process);
          p_signal->SetId();
          p_signal->SetStatus(blob_status::needs_signal);
          m_blobs.push_back(p_signal);
	  weight=1.;
	  pit=p_phases->begin();
	  break;
	default:
	  THROW(fatal_error,"Invalid return value");
	}
	if (weight==0.0 && rpa.gen.NumberOfDicedEvents()==
	    rpa.gen.NumberOfEvents()) return true;
      }
    } while (m_blobs.empty() || p_signal->NInP()==0);
    if (!m_blobs.FourMomentumConservation()) return false;
    {
      Blob *sp(m_blobs.FindFirst(btp::Signal_Process));
      if (sp) {
        double trials((*sp)["Trials"]->Get<double>());
        sp->AddData("Trials",new Blob_Data<double>(trials+m_addn));
        double cxs(m_blobs.Weight());
        m_n+=trials+m_addn;
        m_sum+=cxs;
        m_sumsqr+=sqr(cxs);
        m_addn=0.0;
      }
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) {
	switch ((*pit)->Treat(&m_blobs,weight)) {
	case Return_Value::Nothing :
	  break;
	case Return_Value::Success : 
	  rvalue.IncCall((*pit)->Name());
	  msg_Tracking()<<ATOOLS::om::blue<<"Event_Handler::GenerateEvent("<<mode<<"): "
			<<ATOOLS::om::reset
			<<"Event phase "<<ATOOLS::om::bold<<(*pit)->Name()<<ATOOLS::om::reset
			<<" yields "<<ATOOLS::om::bold<<true<<ATOOLS::om::reset<<std::endl;
	  break;
	case Return_Value::Error :
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncError((*pit)->Name());
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
  case 9000:
  case 9001:
  case 9002: {
    Reset();
    m_blobs.Clear();
    bool flag = false;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::External_MC) {
	bool result=(*pit)->Treat(&m_blobs,weight);
	if (result) flag = true;
      }
    }
    if (!flag) {
      msg_Error()<<"Error in Event_Handler::GenerateEvent"<<std::endl
		 <<"   Treat by External_MC failed, mode = "<<mode<<"."<<std::endl
		 <<"   Return 0 and hope for the best."<<std::endl;
    }
    if (m_blobs.empty()) {
      msg_Out()<<"Potential error in Event_Handler::GenerateEvent"<<std::endl
	       <<"   Empty bloblist would go into analysis !"<<std::endl;
      return false;
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }

    return true;
  }
  case 9999: {
    Reset();
    bool flag(true);
    while (flag) {
      flag = false;
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	if ((*pit)->Type()==eph::Read_In) {
	  bool result=(*pit)->Treat(&m_blobs,weight);
	  if (result) flag = true;
	}
      }
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }
    return true;
  }
  }
  return false;
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
  double xs(TotalXS()), err(TotalErr());
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
