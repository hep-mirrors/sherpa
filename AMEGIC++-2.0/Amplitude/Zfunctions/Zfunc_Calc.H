#ifndef Zfunc_Calc_h
#define Zfunc_Calc_h

#include <vector>
#include "Basic_Func.H"
#include "Lorentz_Function.H"

#define Scalar_Args

namespace AMEGIC {

  class Basic_Sfuncs;

  class zl {
  public:
    enum code {Y      = -1, 
	       Z      = 0, 
	       VVV    = 1,
	       VVS    = 2,
	       SSV    = 3,
	       SSS    = 4,
	       VVVV   = 5,
	       V4     = 6,
	       G4     = 7,
	       V5     = 9,
               V4V3   = 10, 
               VVSS   = 11,
	       VVSS4  = 12,//the generic interaction of 2 vectors and 2 scalars
	       SSSS   = 13,//the generic interaction of 4 scalars
	       AV4    = 20,//Anomalous gauge coupling
	       FFT    = 60,
	       VVT    = 61,
	       SST    = 62,
	       FFVT   = 63,
	       VVVT   = 64,
	       SSST   = 65,
	       FFGS   = 70,
	       VVGS   = 71,
	       SSGS   = 72,
	       FFVGS  = 73,
	       Unknown= 99,
               Triangle = 101,
	       Box      = 102,
	       Pentagon = 103
    }; 
  };
  
  class Zfunc_Calc : public virtual Basic_Func {
  public:
    int                      ncoupl,narg,pn;
    std::vector<Lorentz_Function> lorentzlist;
    zl::code                 type; 

    Zfunc_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {} 
    virtual ~Zfunc_Calc() {}
    virtual Kabbala Do(); 
    
    virtual int GetScalarNumb();
  };
  
  /*!
    \class zl
    \brief Calculator type codes.
  */

  /*!
    \class Zfunc_Calc
    The mother class of all calculator classes.
  */
  /*!
    \fn virtual Kabbala Zfunc_Calc::Do()
    Pure virtual method, must be overloaded!

    Performs the actual calculation in derived classes.
  */
  /*!
    \fn virtual int Zfunc_Calc::GetScalarNumb()
    Virtual method, returns 0.

    Used to return the number of dummy scalar arguments.
  */
  /*!
    \var int Zfunc_Calc::ncoupl
    Number of complex coupling constants
  */
  /*!
    \var int Zfunc_Calc::narg
    Number of arguments
  */
  /*!
    \var int Zfunc_Calc::pn
    Number of internal propagators
  */
  /*!
    \var std::vector<Lorentz_Function> Zfunc_Calc::lorentzlist
    The Lorentz structure of the calculated function
  */
  /*!
    \var zl::code Zfunc_Calc::type
    Unique identifier for the calculator.
  */
}
#endif



























