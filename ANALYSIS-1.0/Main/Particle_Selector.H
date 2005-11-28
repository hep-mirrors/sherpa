#ifndef Particle_Selector_H
#define Particle_Selector_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {  class Particle_Qualifier_Base; }

namespace ANALYSIS {


  class Particle_Selector : public Primitive_Observable_Base {
    std::string  m_inlistname1, m_inlistname2, m_outlistname, m_qualname;
    int m_mode, m_init;
    ATOOLS::Particle_Qualifier_Base * p_qualifier;
  public:
    Particle_Selector(const std::string & inlistname1,
		      const std::string & inlistname2,
		      const std::string & outlistname,int mode);
    Particle_Selector(const std::string &inlist,const std::string &outlist,
		      const std::string &qualname);
    void CreateParticleList();
    void Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout);
    Primitive_Observable_Base * Copy() const;    
    ~Particle_Selector();
  };

} // namespae ANALYSIS
#endif
