#ifndef Splitting_Function_H
#define Splitting_Function_H

#include "Flavour.H"
#include <iostream>

namespace APACIC {
  class Splitting_Function {
  protected:
    const static double s_Nc;
    const static double s_CF;
    const static double s_CA;
    const static double s_TR;
    // Flavours of this particular splitting function
    ATOOLS::Flavour m_flavs[3]; // a -> b c
    // Domain of last intergral
    double m_zmin,m_zmax;
    // Last integral
    double m_lastint;
  public:
    virtual double operator() (double) = 0;
    virtual double GetZ()=0;             

    virtual double GetWeight(double,double,bool)=0;

    virtual double GetCoupling()=0;             
    virtual double GetCoupling(double t)=0;             

    virtual double CrudeInt(double z0);
    virtual double CrudeInt(double,double) = 0;        
    virtual double GetLastInt();

    virtual double Integral(double zmin,double zmax);        

    virtual void Add(Splitting_Function *);
    virtual void SelectOne();
  
    virtual ATOOLS::Flavour & GetFlA();
    virtual ATOOLS::Flavour & GetFlB();
    virtual ATOOLS::Flavour & GetFlC();

    virtual void PrintStat(int mode=0);
  };


  /*!
    \file
    \brief contains the class APACIC::Splitting_Function
  */

  /*! 
    \class Splitting_Function
    \brief is the base class of all splitting functions 

    The class Spltting_Function is a purely virtual class.
    All specific splitting functions are derived from this class, and 
    have to implement the physics.
    This class can be understood as an interface to access an arbitrary 
    splitting function.
  */

  /*!
    \fn virtual double Splitting_Function::operator() (double) 
    \brief returns the exact value of a splitting function at some given z
  */

  /*!
    \fn virtual double Splitting_Function::GetZ()
    \brief returns some z chosen according to logarithmic accuracy 
      dependent on the lower cut off zmin                         
  */

  /*!
    \fn virtual double Splitting_Function::GetCoupling()
    \brief returns the maximal coupling constant
  */

  /*!
    \fn virtual double Splitting_Function::GetCoupling(double t)
    \brief returns the actual coupling constant at scale t
  */

  /*!
    \fn virtual double Splitting_Function::GetWeight(double,double,bool)
    \brief returns the correction weight of the exact splitting function w.r.t. the
      crude (logarithmic) estimate.
  */

  /*!
    \fn virtual double Splitting_Function::CrudeInt(double z0)
    \brief  returns integral of the splitting function to logarithmic
      accuracy with symmetric boundaries zmin=z_0 and zmax=1-z_0  
  */

  /*!
    \fn virtual double Splitting_Function::CrudeInt(double,double)
    \brief returns integral of the splitting function to logarithmic
      accuracy with asymmetric boundaries zmin and zmax               
  */

  /*!
    \fn virtual double Splitting_Function::GetLastInt()
    \brief returns the last value of the integral.
  */

  /*!
    \fn virtual void Splitting_Function::Add(Splitting_Function *)
    \brief adds a splitting function (used in groups only) 
  */

  /*!
    \fn virtual void Splitting_Function::SelectOne()
    \brief select one splitting function (used in Groups only)
  */

  /*!
    \fn virtual ATOOLS::Flavour & Splitting_Function::GetFlA()
    \brief returns flavour \f$a\f$ in the splitting \f$a \to b \; c\f$.
  */

  /*!
    \fn virtual ATOOLS::Flavour & Splitting_Function::GetFlB()
    \brief returns flavour \f$b\f$ in the splitting \f$a \to b \; c\f$.
  */
    
  /*!
    \fn virtual ATOOLS::Flavour & Splitting_Function::GetFlC()
    \brief returns flavour \f$c\f$ in the splitting \f$a \to b \; c\f$.
  */

  /*!
    \fn virtual void PrintStat(int mode=0)
    \brief Prints the actual splitting described. (esp. usefull for Groups)
  */
} 

#endif 
