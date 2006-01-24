//bof
//Version: 3 ADICIC++-0.0/2005/08/18

//ADICIC-specific provider class for global parameters and compile-time
//information.



#ifndef _Information_hpp_
#define _Information_hpp_ _Information_hpp_


#include "Dipole_Flavour.H"





namespace ADICIC {



  class Dipole_Flavour_Info {

  public:

    Dipole_Flavour_Info();
    ~Dipole_Flavour_Info();

    struct {
      const Dipole_Gluon_Base* pkf[22];
      Dipole_Gluon_G g;
    } gluon;

    struct {
      const Dipole_Quark_Base* pkf[7];
      Dipole_Quark_D d;
      Dipole_Quark_U u;
      Dipole_Quark_S s;
      Dipole_Quark_C c;
      Dipole_Quark_B b;
      Dipole_Quark_T t;
    } quark;

    struct {
      const Dipole_Antiquark_Base* pkf[7];
      Dipole_Antiquark_D d;
      Dipole_Antiquark_U u;
      Dipole_Antiquark_S s;
      Dipole_Antiquark_C c;
      Dipole_Antiquark_B b;
      Dipole_Antiquark_T t;
    } antiq;

  };

  extern const Dipole_Flavour_Info info;



}    //eo namespace ADICIC





#endif    //eo _Information_hpp_



//eof
