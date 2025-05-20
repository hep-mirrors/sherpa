#ifndef ALPACA_Tools_Scatter_Merge_Variables_H
#define ALPACA_Tools_Scatter_Merge_Variables_H

#include <iostream>

namespace ALPACA {

  class Scatter_Merge_Variables {
  private :
    double      p_partons;

    //Contains all Variables needed for scatter/merge processes:
    /*
      General:
        -Part in 1
        -Part in 2
        
      Before scatter/merge:
        -Part recoil
        -tau
        -xsec_oe
        -x
        -kT2
        -mg2, Tstar


      After scatter/merge:
        -accepted
        -type: scatter or merge

      If accepted scatter:
        -

      If accepted merge:
        -

    */
    
  public:
    Scatter_Merge_Variables();
    Scatter_Merge_Variables(double partons);
    ~Scatter_Merge_Variables();


  };
}

#endif
