//bof
//Version: 2 ADICIC++-0.0/2004/09/01

//Helper class, carrying the results of the Sudakov calculation.



#ifndef _Sudakov_Result_H_
#define _Sudakov_Result_H_ _Sudakov_Result_H_


#include <enumextra>
#include "Flavour.H"
#include "Sudakov_Strategy.hpp"





namespace ADICIC {



  struct Sudakov_Result {
    Radiation::Group Rad;
    ATOOLS::kf::code Kfc;
    double           P2t;
    double           X1;
    double           X3;
    Sudakov_Result();
    ~Sudakov_Result() {}
  };



}    //eo namespace ADICIC





#endif    //eo _Sudakov_Result_H_



//eof
