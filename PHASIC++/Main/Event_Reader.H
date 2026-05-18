#ifndef PHASIC__Main__Event_Reader_H
#define PHASIC__Main__Event_Reader_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Org/Getter_Function.H"

#include <vector>
#include <string>

namespace PHASIC {

  struct Event_Reader_Key {
    std::vector<std::string> m_files;
  };// end of struct Event_Reader_Key

  class Event_Reader {
  protected:

    ATOOLS::Cluster_Amplitude *p_ampl;
    ATOOLS::NLO_subevt        *p_sub;

    std::vector<std::string> m_files;

    double m_totalxs, m_unitwgt;
    int    m_compute;
    
  public:
    
    typedef ATOOLS::Getter_Function
    <Event_Reader,Event_Reader_Key> Getter_Function;
    
  public:

    inline Event_Reader(const Event_Reader_Key &key):
      m_files(key.m_files), m_compute(1), p_ampl(NULL), p_sub(NULL) {}
    virtual ~Event_Reader() = default;

    virtual void PrintStatistics(std::ostream&) {};

    // Hook for readers that need a one-shot warm-up phase before the
    // first event can be read. Called by Sherpa's cross-section
    // calculation step (i.e. Matrix_Element_Handler::CalculateTotalXSecs)
    // in place of running its own integration. Readers are expected to
    // make this idempotent so that recursive calls from composite
    // processes (e.g. MC@NLO) only trigger the underlying warm-up once.
    virtual void WarmUp() {}

    // Optional teardown hook. Sherpa calls this during shutdown before so
    // readers that own background workers or library resources (e.g. Pepper)
    // can drain in-flight work and release them while their backing libraries
    // (e.g. MPI, Kokkos) are still alive. Must be idempotent: composite
    // processes (e.g. MC@NLO) may forward the call from several places.
    virtual void Finalize() {}

    virtual ATOOLS::Cluster_Amplitude *ReadEvent() = 0;
    virtual void StepBackward() = 0;

    inline double TotalXS() const { return m_totalxs; }
    inline double UnitWeight() const { return m_unitwgt; }

    inline int Compute() const { return m_compute; }

    inline void SetAmpl(ATOOLS::Cluster_Amplitude *ampl) { p_ampl=ampl; }

    inline ATOOLS::NLO_subevt *SubEvt() const { return p_sub; }

  };// end of class Event_Reader
  
}// end of namespace PHASIC

#endif
