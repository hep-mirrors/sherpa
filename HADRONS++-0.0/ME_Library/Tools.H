#ifndef HadTools
#define HadTools
#include "MathTools.H"

namespace HADRONS
{  
  class Tools {
  public:
    static double   Lambda(double,double,double);
    static Complex  BreitWigner( double s, double Mass2, double MassWidth );
    static Complex  BreitWigner( double s, double Mass2, 
				 double Width, double ms );
    static Complex  BreitWigner( double s, double Mass2, 
				 double Width, double ms1, double ms2 );
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms );
    static double   OffShellMassWidth( double s, double Mass2, 
				       double Width, double ms1, double ms2 );
  }; 

}

#endif
