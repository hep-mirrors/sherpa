#ifndef Splitting_Group_H
#define Splitting_Group_H

#include "Splitting_Function.H"
#include "FastList.H"


namespace APACIC {
  // to be replaced with STL vector
  typedef ATOOLS::FastList<Splitting_Function> SplFunList; 
  typedef ATOOLS::FastList<Splitting_Function>::Iterator SplFunIter;

  class Splitting_Group : public Splitting_Function {
  protected:
    SplFunList  m_group;
    Splitting_Function * p_selected;
    // partial sums
    double * p_partsums;
  public:
    /*! 
      Standard Constructor, if a Splitting_Function is specified it will
      be appended to the list of this particular group.
    */
    Splitting_Group(Splitting_Function * spl=0);
    //! Standard Destructor
    virtual ~Splitting_Group();
    //! Summing over all splittings, filling individual values and total sum
    virtual double CrudeInt(double _zmin, double _zmax);
    /*!
      A load of various functions to return values of the particular
      selected function out of the group, see Splitting_Function.H
      for more details.
    */
    double operator()(double z);
    double Integral(double zmin,double zmax);
    double GetZ();
    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool masses);
    ATOOLS::Flavour & GetFlA();
    ATOOLS::Flavour & GetFlB();
    ATOOLS::Flavour & GetFlC();
    //! Adds a splitting function to some exisiting group.
    virtual void Add(Splitting_Function * spl);
    //! Selects one splitting function out of the group.
    virtual void SelectOne();
    void PrintStat(int mode=0);
    
  };

  /*!
    \file
    \brief contains the class APACIC::Splitting_Group
  */

  /*! 
     \class Splitting_Group
     \brief A group of splitting functions 

     If more than one splitting channel is accessible for a flavour.
     all possible splitting fuctions are contained in a single 
     Splitting_Group.

     For instance for a quark a splitting in quark and gluon as well as
     a splitting in quark and photon is possible.
     \f[
       P_q(z) = P_{q \to qg}(z) + P_{q \to q\gamma}(z)
     \f]

     A usual calling sequence is
      -# CrudeInt() -- detemines integrals of the splitting functions 
         at logarithmic accuracy
      -# SelectOne() -- chooses one group member according to the calculated 
                        \em crude integral
      -# GetZ() -- dices \f$z\f$ according to an aproximated version of the 
              selected splitting function with in the domain used for the CrudeInt()
      -# GetWeight(), GetCoupling(), GetFlA(), etc. -- provide all information
          about the selected splitting function and the neccessary rejection weights
          in order to obtain an exact splitting functions
     .
  */

  /*!
   \fn Splitting_Group::Splitting_Group(Splitting_Function * spl=0);
   \brief  Standard Constructor, if a Splitting_Function is specified it will
      be appended to the list of this particular group.
  */

  /*!
    \fn Splitting_Group::~Splitting_Group();
    \brief Destructor
  */

  /*!
    \fn virtual double Splitting_Group::CrudeInt(double _zmin, double _zmax);
    \brief Summing over all splittings, filling individual values and total sum
  */

  /*!
    \fn double Splitting_Group::operator()(double z);
    \brief  A load of various functions to return values of the particular
      selected function out of the group, see Splitting_Function.H
      for more details.
  */

  /*
    double Splitting_Group::GetZ();
    double Splitting_Group::GetCoupling();
    double Splitting_Group::GetCoupling(double t);
    double Splitting_Group::GetWeight(double z,double pt2,bool masses);
    ATOOLS::Flavour & Splitting_Group::GetFlA();
    ATOOLS::Flavour & Splitting_Group::GetFlB();
    ATOOLS::Flavour & Splitting_Group::GetFlC();
  */
  
  /*!
    \fn virtual void Splitting_Group::Add(Splitting_Function * spl);
    \brief Adds a splitting function to the group.
  */

  /*
    \fn virtual void Splitting_Group::SelectOne();
    \brief Selects one splitting function out of the group.
  */

  /*!
    \fn void Splitting_Group::PrintStat(int mode=0);
    \brief prints all splitting functions contained in this 
           group to the screen.
  */
}
#endif



