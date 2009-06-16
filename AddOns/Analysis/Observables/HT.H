#ifndef Analysis_Observables_HT_H
#define Analysis_Observables_HT_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class HT: public Primitive_Observable_Base {  
  private:

    std::string m_reflist;

  public:

    HT(int type,double xmin,double xmax,int nbins,
       const std::string & listname=std::string(""),
       const std::string & reflistname=std::string(""));
    
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class HT

}// end of namespace ANALYSIS

#endif
