#ifndef HT_H
#define HT_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class HT: public Primitive_Observable_Base {  
  public:

    HT(int type,double xmin,double xmax,int nbins,
       const std::string & listname=std::string(""));
    
    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class HT

}// end of namespace ANALYSIS

#endif
