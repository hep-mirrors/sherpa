#ifndef Primitive_Observable_h
#define Primitive_Observable_h

#include "Primitive_Observable_Base.H"
#include "Jet_Finder.H"

namespace ANALYSIS {

  class Multiplicity : public Primitive_Observable_Base {  
  public:
    Multiplicity(int type,double xmin,double xmax,int nbins,
		 const std::string & listname=std::string(""));

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;
  };


}

#endif
