//bof
//Version: 1 ADICIC++-0.0/2004/02/27

//ADICIC-specific provider class for global parameters and compile-time
//information.



#ifndef _Information_hpp_
#define _Information_hpp_ _Information_hpp_


#include "Dipole_Flavour.H"





namespace ADICIC {



  class Info {

  public:

    struct {
      Dipole_Gluon_G g;
    } gluon;

    struct {
      Dipole_Quark_D d;
      Dipole_Quark_U u;
      Dipole_Quark_S s;
      Dipole_Quark_C c;
      Dipole_Quark_B b;
      Dipole_Quark_T t;
    } quark;

    struct {
      Dipole_Antiquark_D d;
      Dipole_Antiquark_U u;
      Dipole_Antiquark_S s;
      Dipole_Antiquark_C c;
      Dipole_Antiquark_B b;
      Dipole_Antiquark_T t;
    } antiquark;

  };

  extern const Info info;



}    //eo namespace ADICIC





#endif    //eo _Information_hpp_



//eof




