#ifndef Normalized_Observable_H
#define Normalized_Observable_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Normalized_Observable: public Primitive_Observable_Base {  
  protected:
    
    ATOOLS::Histogram *p_obs, *p_norm;

  public:

    // constructors
    Normalized_Observable();
    Normalized_Observable(int type,double xmin,double xmax,int nbins,
			      ATOOLS::Selector_Base * sel=NULL);
    Normalized_Observable(const Normalized_Observable & old);

    // destructor
    ~Normalized_Observable();

    // member functions
    Primitive_Observable_Base &
    operator+=(const Primitive_Observable_Base &obs);

    void Reset();
    void Restore(double scale=1.0);
    void EndEvaluation(double scale=1.0);

  };// end of class Normalized_Observable

}// end of namespace ANALYSIS

#endif
