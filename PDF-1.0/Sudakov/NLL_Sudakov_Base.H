#ifndef NLL_Sudakov_Base_H
#define NLL_Sudakov_Base_H

#include "Function_Base.H"
#include "MathTools.H"
#include "Gauss_Integrator.H"
#include "Flavour.H"

namespace ATOOLS {

  class Integration_Info;

}

namespace SHERPA {
  class Sudakov {
  public:
    enum code {
      dontcut     = 0,
      cutatone    = 1,
      cutatkinlim = 2,

      table    = 128,
      numeric  = 256,
      analytic = 384,
      create_table  = 512
    };
  };

  class NLL_Sudakov_Base  : public ATOOLS::Function_Base {

  public:
    virtual double operator()(double Q, double q) = 0;
    virtual double operator()(double) =0;
    virtual double Log(double,double) = 0;
    virtual void AssignKey(ATOOLS::Integration_Info *const info) = 0;
    virtual ~NLL_Sudakov_Base() {};
  };  

  class NLL_Dummy_Sudakov : public NLL_Sudakov_Base {

  public:
    double operator()(double,double) { return 1.; }
    double operator()(double) { return 1.; }
    void   AssignKey(ATOOLS::Integration_Info *info) {};
    double Log(double,double) { return 0; }

  };
  
  /*!
    \file 
    \brief contains the classes SHERPA::NLL_Sudakov_Base and SHERPA::NLL_Dummy_Sudakov
  */
  /*!
    \class NLL_Sudakov_Base
    \brief is a pure virtual base class, providing an interface to any sudakov like object
  */
  /*!
    \class NLL_Dummy_Sudakov
    \brief is a simple example of an sudakov returning always one
           can be used for unknow/not-strong flavours.
  */

}


#endif
