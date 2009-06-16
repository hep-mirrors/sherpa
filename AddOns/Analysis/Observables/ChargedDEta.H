#ifndef Analysis_Observables_ChargedDEta_H
#define Analysis_Observables_ChargedDEta_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class ChargedDEta: public Primitive_Observable_Base {  
    ATOOLS::Flavour m_flav1, m_flav2;
  public:

    ChargedDEta(int type,double xmin,double xmax,int nbins,
                const std::string & listname,
                const int kf1, const int kf2);
    
    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class ChargedDEta

}// end of namespace ANALYSIS

#endif
