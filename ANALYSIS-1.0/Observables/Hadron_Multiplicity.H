#ifndef Hadron_Multiplicity_h
#define Hadron_Multiplicity_h


#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Hadron_Multiplicity : public Primitive_Observable_Base {  
  public:
    Hadron_Multiplicity(const String_Matrix &);
    Hadron_Multiplicity(const Hadron_Multiplicity *);

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;
  };// end of class Hadron_Multiplicity
}// end of namespace ANALYSIS

#endif
