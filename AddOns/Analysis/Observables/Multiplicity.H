#ifndef Analysis_Observables_Multiplicity_H
#define Analysis_Observables_Multiplicity_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "AddOns/Analysis/Observables/Normalized_Observable.H"

namespace ANALYSIS {

  class Multiplicity : public Primitive_Observable_Base {  
  public:
    Multiplicity(int type,double xmin,double xmax,int nbins,
		 const std::string & listname=std::string(""));

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class Multiplicity

  class Hadron_Multiplicities : public Normalized_Observable {  
  public:
    Hadron_Multiplicities(int type,double xmin,double xmax,int nbins,
			  const std::string & listname=std::string(""));

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class Hadron_Multiplicities
}// end of namespace ANALYSIS

#endif
