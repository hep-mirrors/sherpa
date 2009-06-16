#ifndef APACIC_Showers_Splitting_Group_H
#define APACIC_Showers_Splitting_Group_H

#include "APACIC++/Showers/Splitting_Function.H"

namespace APACIC {

  typedef std::vector<double>              Double_Vector; 
  typedef std::vector<Splitting_Function*> Splitting_Vector; 

  class Splitting_Group : public Splitting_Function {
  protected:

    Splitting_Vector m_splittings;
    Double_Vector    m_partsums;

    Splitting_Function *p_selected;

  public:

    // constructor
    Splitting_Group(ATOOLS::Mass_Selector *&ms,
		    Splitting_Function *const spl=NULL);

    // destructor
    virtual ~Splitting_Group();

    // inline functions
    virtual double CrudeInt(double zmin, double zmax);

    double operator()(double z);
    double Integral(double zmin,double zmax);
    double GetZ();
    double GetPhi(double z);

    const ATOOLS::Simple_Polarisation_Info 
    GetPolB(double z, double phi);
    const ATOOLS::Simple_Polarisation_Info 
    GetPolC(double z, double phi, double phib);

    double GetCoupling();
    double GetCoupling(double t);
    double GetWeight(double z,double pt2,bool masses);

    const ATOOLS::Flavour &GetFlA() const;
    const ATOOLS::Flavour &GetFlB() const;
    const ATOOLS::Flavour &GetFlC() const;

    virtual void Add(Splitting_Function * spl);
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



