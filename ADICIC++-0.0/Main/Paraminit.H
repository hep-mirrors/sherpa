//bof
//Version: 2 ADICIC++-0.0/2004/08/10

//The central Dipole_Parameter initialization object.



#ifndef _Paraminit_H_
#define _Paraminit_H_ _Paraminit_H_


//#include "..."


#define PARAMINIT_OUTPUT PARAMINIT_OUTPUT
//#undef  PARAMINIT_OUTPUT





namespace ADICIC {



  class Dipole_Parameter_Init {
  private:
    static unsigned s_status;
  public:
    Dipole_Parameter_Init() {}
    ~Dipole_Parameter_Init() {}
  public:
    static inline const unsigned Status();
    static const bool DoIt();
  };



  inline const unsigned Dipole_Parameter_Init::Status() {    //Static.
    return s_status;
  }



}    //eo namespace ADICIC





#endif    //eo _Paraminit_H_



//eof
