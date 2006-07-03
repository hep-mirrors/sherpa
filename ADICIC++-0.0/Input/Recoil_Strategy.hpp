//bof
//Version: 4 ADICIC++-0.0/2006/05/26

//Possibility of having information at compile time.
//Globally defined parameter sets influencing the calculation of recoils.



#ifndef _Recoil_Strategy_hpp_
#define _Recoil_Strategy_hpp_ _Recoil_Strategy_hpp_


#include <cstdlib>





namespace ADICIC {



  namespace Recoil_Strategy {

    //The settings here must agree with the ones in ADICIC::MakeRecos(..),
    //see Recoil_Calculator.H

    enum Type {
      Unknown    = 0,
      Kleiss     = 2,
      FixDir1    = 1,
      FixDir3    = 3,
      MinimizePt = 4,
      Lonnblad   = 5,
      OldAdicic  = 6,
      Test       = 7,
      Ktii       = 8,
      Kleissii   = 9,
      stop       = -99
    };

    static const size_t NumberOfTypes = 11;

    static const Type List[NumberOfTypes]={
      Unknown,
      Kleiss,
      FixDir1,
      FixDir3,
      MinimizePt,
      Lonnblad,
      OldAdicic,
      Test,
      Ktii,
      Kleissii,
      stop    //Keep this at the bottom.
    };

  }    //eo namespace Recoil_Strategy



}    //eo namespace ADICIC





#endif    //eo _Recoil_Strategy_hpp_



//eof
