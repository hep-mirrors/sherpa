#ifndef AMEGIC_Amplitude_Zfunctions_Basic_Func_H
#define AMEGIC_Amplitude_Zfunctions_Basic_Func_H

#include "ATOOLS/Math/Kabbala.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
 
namespace AMEGIC {

  class Pfunc;
  class Virtual_String_Generator;

  class Spinor {
  public:
    enum code {None = 0,
               u    = 1,
	       ubar = 2,
               v    = 3,
               vbar = 4, 
               Unnown=99
    };
  };

  class Direction {
  public:
    enum code {None     =  0,
	       Incoming = -1,
               Outgoing =  1,
               Unknown  = 99
    };
  }; 

  class Argument {
  public:
    kf_code        kfcode;
    Spinor::code    spinortype;
    Direction::code direction;
    int             numb;
    bool            maped;
    Argument() {
      spinortype = Spinor::None;
      direction  = Direction::Outgoing;
      numb       = -99;
      maped      = false;
    }   
  };

  inline bool operator==(const Argument& a, const Argument& b) 
  {
    if (a.spinortype!=b.spinortype) return false;
    if (a.direction!=b.direction) return false;
    if (a.numb!=b.numb) return false;

    return true;
  }

  inline bool operator!=(const Argument& a, const Argument& b) {return !(a==b);}


  class Basic_Func {
  public:
    int*                      arg;
    Complex*                  coupl;
    Argument*                 ps;
    int                       pn;
    Pfunc_List*               pl;
    Virtual_String_Generator* sgen;
    Basic_Sfuncs*             BS;

    Basic_Func(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : sgen(_sgen), BS(_BS) {}
    virtual ~Basic_Func();

    virtual ATOOLS::Kabbala X(const int,const int,const int);
    virtual ATOOLS::Kabbala V(const int,const int);
    virtual ATOOLS::Kabbala Vcplx(const int,const int,const int s=1);

    void SetArgCouplProp(int,int*,Complex*,int,Argument*,Pfunc_List*);
    void Map(int&);
    void Map(int&,bool&);
    double GetPMass(int,int);
  };

  class Basic_Yfunc : public virtual Basic_Func {
  public:
    Basic_Yfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    

    ATOOLS::Kabbala Y(const int);
    Complex Ycalc(const int,const int,const int,const int,const Complex&,const Complex&);

    template <int,int>
    inline Complex YT(const int,const int,const Complex&,const Complex&);
  };

  class Basic_Zfunc : public virtual Basic_Func {
  public:
    Basic_Zfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    ATOOLS::Kabbala Z(const int,const int);
    Complex Zcalc(const int,const int,const int,const int,
		  const int,const int,const int,const int,
		  const Complex&,const Complex&,const Complex&,const Complex&);

    int Zmassless(const int,const int,const int,const int,
		  const int,const int,const int,const int,
		  const Complex&,const Complex&,const Complex&,const Complex&);
    template <int,int,int,int>
    inline Complex ZT(const int,const int,
		      const int,const int,
		      const Complex&,const Complex&,
		      const Complex&,const Complex&);

    template <int,int,int,int>
    inline Complex ZTM(const int,const int,
		       const int,const int,
		       const Complex&,const Complex&,
		       const Complex&,const Complex&);
  };

  class Basic_Xfunc : public virtual Basic_Func {
  public:
    Basic_Xfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    ATOOLS::Kabbala X(const int,const int);
    ATOOLS::Kabbala X(const int,const int,const int);
    Complex Xcalc(const int,const int,const int,
		  const int,const int,
		  const Complex&,const Complex&);

    template <int,int>
    inline Complex XT(const int,const int,const int,
		      const Complex&,const Complex&);
  };

  class Basic_Mfunc : public virtual Basic_Func {
  public:
    Basic_Mfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    ATOOLS::Kabbala M(const int);
  };

  class Basic_MassTermfunc : public virtual Basic_Func {
  public:
    Basic_MassTermfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    ATOOLS::Kabbala MassTerm(int);
    Complex MassTermCalc(int,int);
    Complex MassTermCalc(int,ATOOLS::Flavour);
  };

  class Basic_Vfunc : public virtual Basic_Func {
  public:
    Basic_Vfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    ATOOLS::Kabbala V(const int,const int);
    Complex Vcalc(const int,const int);
    ATOOLS::Kabbala Vcplx(const int a,const int b,const int s=1);
    Complex Vcplxcalc(const int,const int);
  };

  class Basic_Pfunc : public virtual Basic_Func {
    double Ifunc(double,int);
    double IEfunc(double,int);
    Complex KKProp(double);
  public:
    Basic_Pfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {} 
    Complex Propagator(double,ATOOLS::Flavour);
    Complex Pcalc(const int,const int);
    Complex Pcalc(const ATOOLS::Flavour&,const int);
    ATOOLS::Kabbala P(Pfunc*);
  };

  class Basic_Epsilonfunc : public virtual Basic_Func {
  public:
    double EC(const ATOOLS::Vec4D*,const ATOOLS::Vec4D*,const ATOOLS::Vec4D*,const ATOOLS::Vec4D*);
  public:
    Basic_Epsilonfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : Basic_Func(_sgen,_BS)  {}    
    ATOOLS::Kabbala Epsilon(const int,const int,const int,const int,int=1);
    Complex EpsCalc(const int,const int,const int,const int,const int);
    template <int> Complex EpsCalc(const int,const int,const int,const int);
  };
  template <> Complex Basic_Epsilonfunc::EpsCalc<0>(const int,const int,const int,const int);
  template <> Complex Basic_Epsilonfunc::EpsCalc<1>(const int,const int,const int,const int);
  template <> Complex Basic_Epsilonfunc::EpsCalc<2>(const int,const int,const int,const int);
  template <> Complex Basic_Epsilonfunc::EpsCalc<3>(const int,const int,const int,const int);
  template <> Complex Basic_Epsilonfunc::EpsCalc<4>(const int,const int,const int,const int);

  class Unitarityfunc : public virtual Basic_Func {
    double m_n,m_lambda2;
  public:
    Unitarityfunc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS);
    Complex Ucalc();
    ATOOLS::Kabbala U();
  };


  /*!
    \class Basic_Func
    The mother class of all basic functions and calculator funtions.
  */
  /*!
    \fn void Basic_Func::SetArgCouplProp(int,int*,Complex*,int,Argument*,Pfunc_List*)
    Copies a pointer to arrays of coupling constants, arguments and propagator lists,
    to be used by the derived classes.
  */
  /*!
    \fn void Basic_Func::Map(int&)
    Maps arguments from amplitude labels to the position of the corresponding vector
    in Basic_Sfuncs::Momlist.
  */
  /*!
    \fn void Basic_Func::Map(int&,bool&)
    Maps propagators from amplitude labels to the position of the corresponding vector
    in Basic_Sfuncs::Momlist and ensures no double mapping.
  */
  /*!
    \fn double Basic_Func::GetPMass(int a,int sign)
    Determines the particle mass of propagator with label a, if sign has the mt code for
    a scalar polarization, otherwise 0.
  */


  /*!
    \class Basic_Xfunc
    \brief A basic structure within the helicity formalism.

    The definition for the basic X function is
    \f[
    X(p_1,\lambda_1;p_2;p_3,\lambda_3;c_L,c_R) = 
    \bar u(p_1,\lambda_1) p\!\!\!/_2[c_LP_L+c_RP_R]u(p_3,\lambda_3).
    \f]
    
    Using a spinor basis, the results reads 
    <table border>    
      <tr>
        <td> \f$\lambda_1\lambda_3\f$ </td>
        <td> \f$X(p_1,\lambda_1;p_2;p_3,\lambda_3;c_L,c_R)\f$ </td>
      </tr>
      <tr>
        <td> \f$++\f$ </td>
	<td> \f$ \mu_1\mu_3\eta_2^2c_L+\mu_2^2\eta_1\eta_3c_R 
	+c_RS(+;p_1,p_2)S(-;p_2,p_3)\f$ </td>
      </tr>
      <tr>
        <td> \f$+-\f$ </td>
	<td> \f$c_L\mu_1\eta_2S(+;p_2,p_3)+c_R\mu_3\eta_2S(+;p_1,p_2)\f$ </td>
      </tr>
    </table>

    Missing combinations can be obtained through the 
    replacements \f$+ \leftrightarrow -\f$ and \f$L \leftrightarrow R\f$.
    
    For the definition of \f$\eta\,,\mu\f$ and S functions see class Basic_Sfuncs.
  */
  /*!
    \fn ATOOLS::Kabbala Basic_Xfunc::X(const int a,const int b)
    This method is usually called from calculator classes (Zfunc_Calc).
    a is an index for the argument array Basic_Func::arg, 
    b of the propagator array Basic_Func::ps.

    The method calls the String_Generator and Basic_Xfunc::Xcalc for the numerics.
  */
  /*!
    \fn ATOOLS::Kabbala Basic_Xfunc::X(const int a,const int b,const int m)
    This method is usually called from Basic_Zfuncs, if there is the dummy 99 in the argument 
    list to indicate an explicit polarization vector.

    Calls Basic_Xfunc::Xcalc with the polarization vector in  \f$p_2\f$.
  */
  /*!
    \fn Complex Basic_Xfunc::Xcalc(const int t1,const int s1,const int t2,const int t3,const int s2,const Complex& cR,const Complex& cL)
    Performs the numerical calculation as in the table in the class description.
    t1,t2,t3 are indices for Basic_Sfuncs::Momlist and refer to p; s1,s3 refer to \f$\lambda\f$.
  */  

  /*!
    \class Basic_Yfunc
    \brief A basic structure within the helicity formalism.

    The definition for the basic Y function is
    \f[
    Y(p_1,\lambda_1;p_2,\lambda_2;{c_R,c_L}) \equiv
    \bar u(p_1,\lambda_1)\left[c_R P_R + c_L P_L\right]u(p_2,\lambda_2)
    \f]
    
    Using a spinor basis, the results reads 
    <table border>    
      <tr>
        <td> \f$\lambda_1\lambda_2\f$ </td>
        <td> \f$Y(p_1,\lambda_1;p_2,\lambda_2;c_L,c_R)\f$ </td>
      </tr>
      <tr>
        <td> \f$++\f$ </td>
	<td> \f$ c_R\mu_1\eta_2 + c_L\mu_2\eta_1\f$ </td>
      </tr>
      <tr>
        <td> \f$+-\f$ </td>
	<td> \f$ c_LS(+;p_1,p_2)\f$ </td>
      </tr>
    </table>

    Missing combinations can be obtained through the 
    replacements \f$+ \leftrightarrow -\f$ and \f$L \leftrightarrow R\f$.
    
    For the definition of \f$\eta\,,\mu\f$ and S functions see class Basic_Sfuncs.
  */
  /*!
    \fn ATOOLS::Kabbala Basic_Yfunc::Y(const int a)
    This method is usually called from calculator classes (Zfunc_Calc).
    a is an index of the argument array Basic_Func::arg.

    The method calls the String_Generator and Basic_Yfunc::Ycalc for the numerics.
  */
  /*!
    \fn Complex Basic_Yfunc::Ycalc(const int t1,const int s1,const int t2,const int s2,const Complex& cR,const Complex& cL)
    Performs the numerical calculation as in the table in the class description.
    t1,t2 are indices for Basic_Sfuncs::Momlist and refer to p; s1,s2 refer to \f$\lambda\f$.
  */  

  /*!
    \class Basic_Zfunc
    \brief A basic structure within the helicity formalism.

    The definition for the basic Z function is
    \f[
    Z(p_1,\lambda_1;p_2,\lambda_2; p_3,\lambda_3;p_4,\lambda_4;
    c_L^{12},c_R^{12};c_L^{34},c_R^{34}) =
    \bar u(p_1,\lambda_1)\gamma^\mu [c_L^{12}P_L+c_R^{12}P_R]u(p_2,\lambda_2)
    \bar u(p_3,\lambda_3)\gamma_\mu [c_L^{34}P_L+c_R^{34}P_R]u(p_4,\lambda_4)    
    \f]
    
    Using a spinor basis, the results reads 
    <table border>    
      <tr>
        <td> \f$\lambda_1\lambda_2\lambda_3\lambda_4\f$ </td>
        <td> \f$Z(p_1,\lambda_1;p_2,\lambda_2;
	        p_3,\lambda_3;p_4,\lambda_4; c_L^{12},c_R^{12};c_L^{34},c_R^{34})\f$ </td>
      </tr>
      <tr>
        <td> \f$++++\f$ </td>
	<td> \f$2\left[S(+;p_3,p_1)S(-;p_2,p_4)c_R^{12}c_R^{34}+ \mu_1\mu_2\eta_3
	        \eta_4c_L^{12}c_R^{34}+ \mu_3\mu_4\eta_1\eta_2c_R^{12}c_L^{34}\right]\f$ </td>
      </tr>
      <tr>
        <td> \f$+++-\f$ </td>
	<td> \f$2\eta_2c_R^{12}\left[S(+;p_1,p_4)\mu_3c_L^{34} +S(+;p_1,p_3)
	     \mu_4c_R^{34}\right]\f$ </td>
      </tr>
      <tr>
        <td> \f$++-+\f$ </td>
	<td> \f$2\eta_1c_R^{12}\left[S(-;p_3,p_2)\mu_4c_L^{34} q +S(-;p_4,p_2)
	     \mu_3c_R^{34}\right]\f$ </td>
      </tr>
      <tr>
        <td> \f$++--\f$ </td>
	<td> \f$2\left[S(+;p_4,p_1)S(-;p_2,p_3)c_R^{12}c_L^{34}+ \mu_1\mu_2\eta_3
	        \eta_4c_L^{12}c_L^{34}+ \mu_3\mu_4\eta_1\eta_2c_R^{12}c_R^{34}\right]\f$ </td>
      </tr>
      <tr>
        <td> \f$+-++\f$ </td>
	<td> \f$2\eta_4c_R^{34}\left[S(+;p_1,p_3)\mu_2c_R^{12} +S(+;p_2,p_3)
	        \mu_1c_L^{12}\right]\f$ </td>
      </tr>
      <tr>
        <td> \f$+-+-\f$ </td>
	<td> 0 </td>
      </tr>
      <tr>
        <td> \f$+--+\f$ </td>
	<td> \f$-2\left[ \mu_1\mu_4\eta_2\eta_3c_L^{12}c_L^{34} +\mu_2\mu_3
	        \eta_1\eta_4c_R^{12}c_R^{34} -\mu_1\mu_3\eta_2\eta_4c_L^{12}c_R^{34}
	        -\mu_2\mu_4\eta_1\eta_3c_R^{12}c_L^{34}\right]\f$ </td>
      </tr>
      <tr>
        <td> \f$+---\f$ </td>
	<td> \f$2\eta_3c_L^{34}\left[S(+;p_4,p_2)\mu_1c_L^{12} +S(+;p_1,p_4)
	        \mu_2c_R^{12}\right]\f$ </td>
      </tr>
    </table>

    Missing combinations can be obtained through the 
    replacements \f$+ \leftrightarrow -\f$ and \f$L \leftrightarrow R\f$.
    
    For the definition of \f$\eta\,,\mu\f$ and S functions see class Basic_Sfuncs.
  */
  /*!
    \fn ATOOLS::Kabbala Basic_Zfunc::Z(const int a, const int b)
    This method is usually called from calculator classes (Zfunc_Calc).
    a,b are indices of the argument array Basic_Func::arg.

    The method calls the String_Generator and Basic_Zfunc::Zcalc for the numerics.
  */
  /*!
    \fn Complex Basic_Zfunc::Zcalc(const int t1,const int s1,const int t2,const int s2,const int t3,const int s3,const int t4,const int s4,const Complex& cR1,const Complex& cL1,const Complex& cR2,const Complex& cL2)
    Performs the numerical calculation as in the table in the class description.
    t1-t4 are indices for Basic_Sfuncs::Momlist and refer to p; s1-s4 refer to \f$\lambda\f$.
  */  

  /*!
    \class Basic_Vfunc
    A basic function to calculate a scalar product of four-vectors 
  */
  /*!
    \fn ATOOLS::Kabbala Basic_Vfunc::V(const int a,const int b)
    This method is usually called from calculator classes (Zfunc_Calc).
    a,b are indices of the propagator array Basic_Func::ps.

    The method calls the String_Generator and Basic_Vfunc::Vcalc for the numerics.
  */
  /*!
    \fn Complex Basic_Vfunc::Vcalc(const int a,const int b)
    Numerical calculation of the scalar product.
    a,b are indices for Basic_Sfuncs::Momlist.
  */  
  /*!
    \fn ATOOLS::Kabbala Basic_Vfunc::Vcplx(const int a,const int b,const int s)
    Polarization vectors may be complex. This method checks the type of a vector
    and calculates it using Basic_Vfunc::Vcalc or Basic_Vfunc::Vcplxcalc respectively.
    a,b are indices for Basic_Sfuncs::Momlist.

    The method calls the String_Generator and Basic_Vfunc::Vcalc for the numerics.
  */
  /*!
    \fn Complex Basic_Vfunc::Vcplxcalc(const int a,const int b)
    Numerical calculation for a complex scalar product.
    a,b are indices for Basic_Sfuncs::Momlist.
  */  

  /*!
    \class Basic_Mfunc
    Basic function to determine the mass term in massive propagators.
  */
  /*!
    \fn ATOOLS::Kabbala Basic_Mfunc::M(const int a)
    Returns \f$\frac{1}{\rm{Mass}^2}\f$ for the propagator with index a in Basic_Func::ps.
    If the propagator is massless 0 is returned.
  */

  /*!
    \class Basic_Pfunc
    Basic function to calculate the propagator factors. 
  */
}

#include "Basic_Func.icc"
//#include "AMEGIC++/String/String_Generator.H"


#endif
