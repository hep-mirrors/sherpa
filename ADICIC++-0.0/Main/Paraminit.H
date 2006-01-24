//bof
//Version: 3 ADICIC++-0.0/2005/06/07

//The central Dipole_Parameter initialization object.



#ifndef _Paraminit_H_
#define _Paraminit_H_ _Paraminit_H_


//#include "..."


#define PARAMINIT_OUTPUT PARAMINIT_OUTPUT
//#undef  PARAMINIT_OUTPUT





namespace ADICIC {



  class Dipole_Parameter_Init {
  private:
    static const bool Do();
  public:
    Dipole_Parameter_Init();
    ~Dipole_Parameter_Init() {}
  };



}    //eo namespace ADICIC





#endif    //eo _Paraminit_H_



//eof
