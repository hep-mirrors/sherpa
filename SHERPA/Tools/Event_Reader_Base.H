#ifndef SHERPA_Tools_Event_Reader_Base_H
#define SHERPA_Tools_Event_Reader_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle.H"


namespace SHERPA {
  class Event_Reader_Base {
  protected:
    std::string        m_path, m_add, m_file;

    double             m_weight;

  public:
    Event_Reader_Base(const std::string & path,const std::string & file)
      : m_path(path), m_add(file), m_file(file) {}
    virtual ~Event_Reader_Base() {}
    virtual bool FillBlobs(ATOOLS::Blob_List *)=0;

    double      Weight()  { return m_weight; }
    std::string GetPath() { return m_path; }
  };
}

#endif
