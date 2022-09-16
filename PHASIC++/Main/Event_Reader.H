#ifndef PHASIC__Main__Event_Reader_H
#define PHASIC__Main__Event_Reader_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Getter_Function.H"

#include <vector>
#include <string>

namespace PHASIC {

  struct Event_Reader_Key {
    std::vector<std::string> m_files;
  };// end of struct Event_Reader_Key

  class Event_Reader {
  protected:

    std::vector<std::string> m_files;

    double m_totalxs;
    bool   m_compute;
    
  public:
    
    typedef ATOOLS::Getter_Function
    <Event_Reader,Event_Reader_Key> Getter_Function;
    
  public:

    inline Event_Reader(const Event_Reader_Key &key):
      m_files(key.m_files), m_compute(true) {}
    virtual ~Event_Reader() = default;

    virtual ATOOLS::Cluster_Amplitude *ReadEvent() = 0;

    inline double TotalXS() const { return m_totalxs; }

    inline bool Compute() const { return m_compute; }

  };// end of class Event_Reader
  
}// end of namespace PHASIC

#endif
