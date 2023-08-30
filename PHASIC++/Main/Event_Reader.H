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

    virtual ATOOLS::Cluster_Amplitude *ReadEvent() = 0;

    inline double TotalXS() const { return m_totalxs; }
    inline double UnitWeight() const { return m_unitwgt; }

    inline int Compute() const { return m_compute; }

    inline void SetAmpl(ATOOLS::Cluster_Amplitude *ampl) { p_ampl=ampl; }

    inline ATOOLS::NLO_subevt *SubEvt() const { return p_sub; }

  };// end of class Event_Reader
  
}// end of namespace PHASIC

#endif
