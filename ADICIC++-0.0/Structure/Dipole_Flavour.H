//bof
//Version: 1 ADICIC++-0.0/2004/06/05

//Purpose is to have defined dipole-flavour structures throughout the model.



#ifndef _Dipole_Flavour_H_
#define _Dipole_Flavour_H_ _Dipole_Flavour_H_


#include <enumextra>
#include "Flavour.H"





namespace ADICIC {



  class Dipole_Flavour_Init {
  private:
    static const bool status;
    static const bool Now() {
      ATOOLS::ParticleInit("./");
    }
  public:
    Dipole_Flavour_Init() {}
    ~Dipole_Flavour_Init() {}
  };










  class Dipole_Gluon_Base {
  protected:
    static const short GNtag;
    static const ATOOLS::Flavour GN;

    Dipole_Gluon_Base(const Dipole_Gluon_Base&) {}
    Dipole_Gluon_Base& operator=(const Dipole_Gluon_Base&) { return *this;}

  public:
    Dipole_Gluon_Base() {}
    virtual ~Dipole_Gluon_Base() {}

    virtual const short Tag() const { return GNtag;}
    virtual const ATOOLS::Flavour operator()() const=0;
    //virtual const ATOOLS::Flavour& operator()() const=0;
  };





  class Dipole_Gluon_G : public Dipole_Gluon_Base {
  private:
    Dipole_Gluon_G(const Dipole_Gluon_G&) {}
    Dipole_Gluon_G& operator=(const Dipole_Gluon_G&) { return *this;}
  public:
    Dipole_Gluon_G();
    ~Dipole_Gluon_G();

    //If several gluon-type-like particles exist then overwriting is necessary.
    //const short Tag() const { return GNtag;}    //Special case!

    const ATOOLS::Flavour operator()() const { return GN;}
  };










  class Dipole_Quark_Base {
  protected:
    static const short QUtag;
    static const ATOOLS::Flavour DQ;
    static const ATOOLS::Flavour UQ;
    static const ATOOLS::Flavour SQ;
    static const ATOOLS::Flavour CQ;
    static const ATOOLS::Flavour BQ;
    static const ATOOLS::Flavour TQ;

    Dipole_Quark_Base(const Dipole_Quark_Base&) {}
    Dipole_Quark_Base& operator=(const Dipole_Quark_Base&) { return *this;}

  public:
    Dipole_Quark_Base() {}
    virtual ~Dipole_Quark_Base() {}

    virtual const short Tag() const { return QUtag;}
    virtual const ATOOLS::Flavour operator()() const=0;
    //virtual const ATOOLS::Flavour& operator()() const=0;
  };





  class Dipole_Quark_D : public Dipole_Quark_Base {
  private:
    Dipole_Quark_D(const Dipole_Quark_D&) {}
    Dipole_Quark_D& operator=(const Dipole_Quark_D&) { return *this;}
  public:
    Dipole_Quark_D();
    ~Dipole_Quark_D();
    const ATOOLS::Flavour operator()() const { return DQ;}
  };





  class Dipole_Quark_U : public Dipole_Quark_Base {
  private:
    Dipole_Quark_U(const Dipole_Quark_U&) {}
    Dipole_Quark_U& operator=(const Dipole_Quark_U&) { return *this;}
  public:
    Dipole_Quark_U();
    ~Dipole_Quark_U();
    const ATOOLS::Flavour operator()() const { return UQ;}
  };





  class Dipole_Quark_S : public Dipole_Quark_Base {
  private:
    Dipole_Quark_S(const Dipole_Quark_S&) {}
    Dipole_Quark_S& operator=(const Dipole_Quark_S&) { return *this;}
  public:
    Dipole_Quark_S();
    ~Dipole_Quark_S();
    const ATOOLS::Flavour operator()() const { return SQ;}
  };





  class Dipole_Quark_C : public Dipole_Quark_Base {
  private:
    Dipole_Quark_C(const Dipole_Quark_C&) {}
    Dipole_Quark_C& operator=(const Dipole_Quark_C&) { return *this;}
  public:
    Dipole_Quark_C();
    ~Dipole_Quark_C();
    const ATOOLS::Flavour operator()() const { return CQ;}
  };





  class Dipole_Quark_B : public Dipole_Quark_Base {
  private:
    Dipole_Quark_B(const Dipole_Quark_B&) {}
    Dipole_Quark_B& operator=(const Dipole_Quark_B&) { return *this;}
  public:
    Dipole_Quark_B();
    ~Dipole_Quark_B();
    const ATOOLS::Flavour operator()() const { return BQ;}
  };





  class Dipole_Quark_T : public Dipole_Quark_Base {
  private:
    Dipole_Quark_T(const Dipole_Quark_T&) {}
    Dipole_Quark_T& operator=(const Dipole_Quark_T&) { return *this;}
  public:
    Dipole_Quark_T();
    ~Dipole_Quark_T();
    const ATOOLS::Flavour operator()() const { return TQ;}
  };










  class Dipole_Antiquark_Base {
  protected:
    static const short AQtag;
    static const ATOOLS::Flavour DA;
    static const ATOOLS::Flavour UA;
    static const ATOOLS::Flavour SA;
    static const ATOOLS::Flavour CA;
    static const ATOOLS::Flavour BA;
    static const ATOOLS::Flavour TA;

    Dipole_Antiquark_Base(const Dipole_Antiquark_Base&) {}
    Dipole_Antiquark_Base& operator=(const Dipole_Antiquark_Base&) {
      return *this;}

  public:
    Dipole_Antiquark_Base() {}
    virtual ~Dipole_Antiquark_Base() {}

    virtual const short Tag() const { return AQtag;}
    virtual const ATOOLS::Flavour operator()() const=0;
    //virtual const ATOOLS::Flavour& operator()() const=0;
  };
 




  class Dipole_Antiquark_D : public Dipole_Antiquark_Base {
  private:
    Dipole_Antiquark_D(const Dipole_Antiquark_D&) {}
    Dipole_Antiquark_D& operator=(const Dipole_Antiquark_D&) { return *this;}
  public:
    Dipole_Antiquark_D();
    ~Dipole_Antiquark_D();
    const ATOOLS::Flavour operator()() const { return DA;}
  };





  class Dipole_Antiquark_U : public Dipole_Antiquark_Base {
  private:
    Dipole_Antiquark_U(const Dipole_Antiquark_U&) {}
    Dipole_Antiquark_U& operator=(const Dipole_Antiquark_U&) { return *this;}
  public:
    Dipole_Antiquark_U();
    ~Dipole_Antiquark_U();
    const ATOOLS::Flavour operator()() const { return UA;}
  };





  class Dipole_Antiquark_S : public Dipole_Antiquark_Base {
  private:
    Dipole_Antiquark_S(const Dipole_Antiquark_S&) {}
    Dipole_Antiquark_S& operator=(const Dipole_Antiquark_S&) { return *this;}
  public:
    Dipole_Antiquark_S();
    ~Dipole_Antiquark_S();
    const ATOOLS::Flavour operator()() const { return SA;}
  };





  class Dipole_Antiquark_C : public Dipole_Antiquark_Base {
  private:
    Dipole_Antiquark_C(const Dipole_Antiquark_C&) {}
    Dipole_Antiquark_C& operator=(const Dipole_Antiquark_C&) { return *this;}
  public:
    Dipole_Antiquark_C();
    ~Dipole_Antiquark_C();
    const ATOOLS::Flavour operator()() const { return CA;}
  };





  class Dipole_Antiquark_B : public Dipole_Antiquark_Base {
  private:
    Dipole_Antiquark_B(const Dipole_Antiquark_B&) {}
    Dipole_Antiquark_B& operator=(const Dipole_Antiquark_B&) { return *this;}
  public:
    Dipole_Antiquark_B();
    ~Dipole_Antiquark_B();
    const ATOOLS::Flavour operator()() const { return BA;}
  };





  class Dipole_Antiquark_T : public Dipole_Antiquark_Base {
  private:
    Dipole_Antiquark_T(const Dipole_Antiquark_T&) {}
    Dipole_Antiquark_T& operator=(const Dipole_Antiquark_T&) { return *this;}
  public:
    Dipole_Antiquark_T();
    ~Dipole_Antiquark_T();
    const ATOOLS::Flavour operator()() const { return TA;}
  };



}    //eo namespace ADICIC





#endif    //eo _Dipole_Flavour_H_



//eof







