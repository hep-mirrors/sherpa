#ifndef Particle_Selector_H
#define Particle_Selector_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {


  class Particle_Selector : public Primitive_Observable_Base {
    std::string  m_inlistname, m_outlistname;
    int m_mode;
    ATOOLS::Particle_Qualifier_Base * p_qualifier;
  public:
    Particle_Selector(const std::string & inlistname,const std::string & outlistname, int mode=0);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout);
    Primitive_Observable_Base * Copy() const;    
    ~Particle_Selector();
  };

} // namespae ANALYSIS
#endif
