#ifndef Universal_Selector_H
#define Universal_Selector_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Universal_Selector: public Primitive_Observable_Base {  
  private:
    std::string m_outlist;
    std::string m_key;
    double m_keymin, m_keymax;
  public:
    Universal_Selector(const std::string & obskey, const double keymin, const double keymax, 
		       const std::string & inlistname, const std::string & outlistname);
    void Evaluate(const ATOOLS::Blob_List & ,double weight, int ncout);
    void CreateParticleList(bool force=true);
    Primitive_Observable_Base *Copy() const;

  };// end of class Universal_Selector

  class General_Observable : public Primitive_Observable_Base {  
    std::string m_key;
  public:

    General_Observable(int type,double xmin,double xmax,int nbins,
		 const std::string & key);

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class General_Observable

}// end of namespace ANALYSIS

#endif
