//bof
//Version: 2 ADICIC++-0.0/2004/09/01

//ADICIC-specific provider class for global parameters and compile-time
//information.



#ifndef _Information_hpp_
#define _Information_hpp_ _Information_hpp_


#include <map>
#include "Dipole_Flavour.H"





namespace ADICIC {



  class Dipole_Flavour_Info {

  public:

    Dipole_Flavour_Info();
    ~Dipole_Flavour_Info();

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
    } antiq;

  };

  extern const Dipole_Flavour_Info info;










  class Dipquarkbox {
  private:
    Dipquarkbox(const Dipquarkbox&);
    Dipquarkbox& operator=(const Dipquarkbox&);
    std::map<ATOOLS::kf::code,const Dipole_Quark_Base*> m_kfqkbase;
  public:
    Dipquarkbox();
    ~Dipquarkbox() {}
    const Dipole_Quark_Base& operator[](ATOOLS::kf::code kfc) {
      return *(m_kfqkbase[kfc]);}
  };
  class Dipantiqbox {
  private:
    Dipantiqbox(const Dipantiqbox&);
    Dipantiqbox& operator=(const Dipantiqbox&);
    std::map<ATOOLS::kf::code,const Dipole_Antiquark_Base*> m_kfaqbase;
  public:
    Dipantiqbox();
    ~Dipantiqbox() {}
    const Dipole_Antiquark_Base& operator[](ATOOLS::kf::code kfc) {
      return *(m_kfaqbase[kfc]);}
  };





  class Dipole_Flavour_Interface {
  private:
    Dipquarkbox m_qk;
    Dipantiqbox m_aq;
  public:
    Dipole_Flavour_Interface();
    ~Dipole_Flavour_Interface();
    Dipquarkbox& quark;
    Dipantiqbox& antiq;
  };

  extern const Dipole_Flavour_Interface interface;



}    //eo namespace ADICIC





#endif    //eo _Information_hpp_



//eof
