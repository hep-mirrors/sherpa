#ifndef Multiplicity_h
#define Multiplicity_h


#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Multiplicity : public Primitive_Observable_Base {  
  public:
    Multiplicity(const String_Matrix &);
    Multiplicity(const Multiplicity *);

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;
  };// end of class Multiplicity
}// end of namespace ANALYSIS

#endif
