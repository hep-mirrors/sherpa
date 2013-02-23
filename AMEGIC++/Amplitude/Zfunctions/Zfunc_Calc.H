#ifndef AMEGIC_Amplitude_Zfunctions_Zfunc_Calc_H
#define AMEGIC_Amplitude_Zfunctions_Zfunc_Calc_H

#include <vector>
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "MODEL/Interaction_Models/Lorentz_Function.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Smart_Pointer.H"

#define Scalar_Args

namespace AMEGIC {

  class Basic_Sfuncs;
  class Zfunc_Generator;
  class Point;
  class Zfunc;

  struct ZFCalc_Key {
    Virtual_String_Generator *p_sgen;
    Basic_Sfuncs *p_bs;
    MODEL::Interaction_Model_Base *p_model;
    ZFCalc_Key(Virtual_String_Generator *const sgen,Basic_Sfuncs *const bs,
	       MODEL::Interaction_Model_Base *const model):
      p_sgen(sgen), p_bs(bs), p_model(model) {}
  };

  class Zfunc_Calc : public virtual Basic_Func {
  public:
    int                      ncoupl,narg,pn;
    std::vector<MODEL::Lorentz_Function*> lorentzlist;
    std::string type; 

    Zfunc_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {} 
    virtual ~Zfunc_Calc();

    Zfunc_Calc *GetCopy() const;
    virtual ATOOLS::Kabbala Do(); 
    
    virtual int GetScalarNumb();
    virtual void SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			 Point *const p,Point *const pf,Point *&pb,
			 int *lfnumb,int *canumb);
  };
  
  typedef ATOOLS::Getter_Function<Zfunc_Calc,ZFCalc_Key> ZFCalc_Getter;

  /*!
    \class zl
    \brief Calculator type codes.
  */

  /*!
    \class Zfunc_Calc
    The mother class of all calculator classes.
  */
  /*!
    \fn virtual ATOOLS::Kabbala Zfunc_Calc::Do()
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

#define DEFINE_ZFCALC_GETTER(CLASS,TAG,INFO)				\
  DECLARE_GETTER(CLASS,TAG,Zfunc_Calc,ZFCalc_Key);			\
  Zfunc_Calc *ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::		\
  operator()(const ZFCalc_Key &key) const				\
  { return new CLASS(key.p_sgen,key.p_bs); }				\
  void ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::			\
  PrintInfo(std::ostream &str,const size_t width) const			\
  { str<<INFO; }
#define DEFINE_ZFLOOPCALC_GETTER(CLASS,TAG,INFO)			\
  DECLARE_GETTER(CLASS,TAG,Zfunc_Calc,ZFCalc_Key);			\
  Zfunc_Calc *ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::		\
  operator()(const ZFCalc_Key &key) const				\
  { if (key.p_model && !key.p_model->HasLoops()) return NULL;		\
    return new CLASS(key.p_sgen,key.p_bs); }				\
  void ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::			\
  PrintInfo(std::ostream &str,const size_t width) const			\
  { str<<INFO; }
#define DEFINE_ZFTENSORCALC_GETTER(CLASS,TAG,INFO)			\
  DECLARE_GETTER(CLASS,TAG,Zfunc_Calc,ZFCalc_Key);			\
  Zfunc_Calc *ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::		\
  operator()(const ZFCalc_Key &key) const				\
  { if (key.p_model && !key.p_model->HasTensors()) return NULL;		\
    return new CLASS(key.p_sgen,key.p_bs); }				\
  void ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::			\
  PrintInfo(std::ostream &str,const size_t width) const			\
  { str<<INFO; }
#define DEFINE_ZFAGCCALC_GETTER(CLASS,TAG,INFO)				\
  DECLARE_GETTER(CLASS,TAG,Zfunc_Calc,ZFCalc_Key);			\
  Zfunc_Calc *ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::		\
  operator()(const ZFCalc_Key &key) const				\
  { if (key.p_model && !key.p_model->HasAGCs()) return NULL;		\
    return new CLASS(key.p_sgen,key.p_bs); }				\
  void ATOOLS::Getter<Zfunc_Calc,ZFCalc_Key,CLASS>::			\
  PrintInfo(std::ostream &str,const size_t width) const			\
  { str<<INFO; }

#endif



























