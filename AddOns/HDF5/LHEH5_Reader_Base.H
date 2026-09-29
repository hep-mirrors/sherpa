#ifndef AddOns__HDF5__LHEH5_Reader_Base_H
#define AddOns__HDF5__LHEH5_Reader_Base_H

// Behaviour shared by every reader of the LHEH5 parton-level event format.
// The concrete readers differ only in where a buffer of events comes from:
// a file on disk for HDF5_Reader, an in-memory database filled by the
// Pepper library for Pepper_Reader. Everything downstream of "here is a
// buffer of events" -- translating one event into a Cluster_Amplitude,
// serving it as often as its trial count demands, and moving on to the
// next buffer once the current one is drained -- lives here.

#include <string>

#include "AddOns/HDF5/LHEH5.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "PHASIC++/Main/Event_Reader.H"

namespace LHEH5 {

  class LHEH5_Reader_Base: public PHASIC::Event_Reader,
                           public ATOOLS::MPI_Object {
  protected:

    // Index of the next event to serve from the current buffer, and the
    // number of times the event held in p_ampl still has to be handed out
    // (an event carries the number of trials that produced it).
    size_t m_ievt, m_trials;
    // Number of events a single buffer holds at most.
    size_t m_ncache;

    ATOOLS::Vec4D_Vector m_ctmoms;

    // The buffer that events are currently served from. Only called while
    // m_ievt < NEventsInBuffer(), i.e. never on an empty buffer.
    virtual const LHEFile &CurrentFile() const = 0;
    virtual size_t NEventsInBuffer() const = 0;

    // Make the next buffer of events the current one, resetting m_ievt.
    // Returns false once no further events can be obtained, which ends
    // event generation with ExhaustionMessage() shown to the user.
    virtual bool NextBuffer() = 0;
    virtual std::string ExhaustionMessage() const = 0;

    bool FillAmplitude()
    {
      const LHEFile &file(CurrentFile());
      Event e(file.GetEvent(m_ievt));
      if (file.Version()[0]==2 &&
	  file.Version()[1]==0 &&
	  file.Version()[2]==0) e.pinfo.unitwgt*=ATOOLS::rpa->Picobarn();
      msg_Debugging()<<e<<"\n";
      if (e.empty()) {
	m_trials+=e.trials;
	m_ievt++;
	return false;
      }
      if (e[0].pz<0 && e[1].pz>0) {
	msg_Debugging()<<"Flip initial states\n";
	std::swap<Particle>(e[0],e[1]);
	std::swap<double>(e.z1,e.z2);
	if (e.ctparts.size()) {
	  std::swap<Particle>(e.ctparts[0],e.ctparts[1]);
	  if (e.ijt<2) e.ijt=1-e.ijt;
	  if (e.kt<2) e.kt=1-e.kt;
	  if (e.i<2) e.i=1-e.i;
	  if (e.k<2) e.k=1-e.k;
	}
      }
      p_ampl = ATOOLS::Cluster_Amplitude::New();
      for (size_t i(0);i<e.size();++i) {
	ATOOLS::Flavour fl((long int)(e[i].id));
	ATOOLS::Vec4D p(e[i].e,e[i].px,e[i].py,e[i].pz);
	ATOOLS::ColorID cl(i<2?e[i].cl2:e[i].cl1,
			   i<2?e[i].cl1:e[i].cl2);
	p_ampl->CreateLeg(i<2?-p:p,i<2?fl.Bar():fl,cl);
      }
      p_ampl->SetNIn(2);
      p_ampl->SetMuR2(ATOOLS::sqr(e.mur));
      p_ampl->SetMuF2(ATOOLS::sqr(e.muf));
      p_ampl->SetMuQ2(ATOOLS::sqr(e.muq));
      p_ampl->SetKT2(ATOOLS::sqr(e.muq));
      p_ampl->SetLKF(e.wgts[0]);
      m_compute=1;
      m_trials+=e.trials;
      if (e.pinfo.npnlo>0) {
	p_ampl->SetLKF(e.psw);
	m_compute=2;
      }
      if (e.ctparts.empty()) p_sub->m_n=0;
      else {
	m_ctmoms.resize(p_sub->m_n=e.ctparts.size());
	for (size_t i(0);i<e.ctparts.size();++i)
	  m_ctmoms[i]=ATOOLS::Vec4D(e.ctparts[i].e,e.ctparts[i].px,
				    e.ctparts[i].py,e.ctparts[i].pz);
	p_sub->p_mom=&m_ctmoms[0];
	p_sub->m_ijt=e.ijt;
	p_sub->m_kt=e.kt;
	p_sub->m_i=e.i;
	p_sub->m_j=e.j;
	p_sub->m_k=e.k;
	p_sub->m_x1=e.z1;
	p_sub->m_x2=e.z2;
	p_sub->m_result=e.bbpsw;
      }
      return true;
    }

  public:

    // `cachesizekey` names the run-card setting that sizes one buffer,
    // `sourcedesc` describes the event source in user-facing messages.
    LHEH5_Reader_Base(const PHASIC::Event_Reader_Key &key,
		      const std::string &cachesizekey,
		      const std::string &sourcedesc):
      Event_Reader(key), m_ievt(0), m_trials(0), m_ncache(0)
    {
      ATOOLS::Settings& s {ATOOLS::Settings::GetMainSettings()};
      m_ncache = s[cachesizekey].SetDefault(10000).Get<int>();

      p_sub = new ATOOLS::NLO_subevt();
      s_objects.push_back(this);

      // Guard against a setting that will, in general, cause synchronization
      // to break in MPI runs, such that the entire run stalls.
      if (s["PRINT_MPI_XS"].Get<bool>() && ATOOLS::mpi->MySize() > 1) {
        THROW(invalid_input,
              "`PRINT_MPI_XS: true` can not be used when reading parton-level "
              "events from " + sourcedesc
              + ". Please use `PRINT_MPI_XS: false`.");
      }
    }

    ~LHEH5_Reader_Base()
    {
      if (p_ampl) p_ampl->Delete();
      delete p_sub;
    }

    void StepBackward() override final
    {
      --m_ievt;
      if (!FillAmplitude())
	THROW(fatal_error, "Could not fill amplitude.");
      m_trials = 1;
    }

    ATOOLS::Cluster_Amplitude *ReadEvent() override final
    {
      DEBUG_FUNC("i="<<m_ievt<<"("<<NEventsInBuffer()<<"),trial="<<m_trials);
      if (m_trials==0 && p_ampl!=NULL) {
	++m_trials;
	--m_ievt;
      }
      if (m_ievt>=NEventsInBuffer()) {
	if (!NextBuffer()) {
	  msg_Info()<<ATOOLS::om::brown<<ATOOLS::om::bold<<"WARNING"
		    <<ATOOLS::om::reset<<": "<<ExhaustionMessage()
		    <<" Will stop event generation ...\n";
	  ATOOLS::rpa->gen.SetNumberOfEvents(
	      ATOOLS::rpa->gen.NumberOfGeneratedEvents());
	  return NULL;
	}
      }
      if (p_ampl==NULL) {
	if (!FillAmplitude())
	  return NULL;
      }
      if (m_trials==1) {
	ATOOLS::Cluster_Amplitude *ampl(p_ampl);
	p_ampl=NULL;
	--m_trials;
	m_ievt++;
	return ampl;
      }
      --m_trials;
      return NULL;
    }

  };// end of class LHEH5_Reader_Base

}// end of namespace LHEH5

#endif
