#ifndef AMISIC_Model_Profile_Function_H
#define AMISIC_Model_Profile_Function_H
/*!
  \file Profile_Function.H
  \brief Declares the classes Flat_Profile, Exponential_Profile,
  Gaussian_Profile and Double_Gaussian_Profile    
*/

#include "AMISIC++/Model/Profile_Function_Base.H"

namespace AMISIC {

  class Flat_Profile: public Profile_Function_Base {
  private:

    double m_radius;

  public:
  
    // constructor
    Flat_Profile(const double radius);

    // member functions
    virtual double Value(const double b) const;
    virtual double MajorValue(const double b) const;

    virtual double MajorIntegral(const double b) const;
    virtual double InverseMajorIntegral(const double I) const;

  };// end of class Flat_Profile
  /*!
    \class Flat_Profile
    \brief Represents a flat matter distribution of the hadron 
    in impact parameter space
  */

  class Exponential_Profile: public Profile_Function_Base {
  private:

    double m_radius;

  public:
  
    // constructor
    Exponential_Profile(const double radius);

    // member functions
    virtual double Value(const double b) const;
    virtual double MajorValue(const double b) const;

    virtual double MajorIntegral(const double b) const;
    virtual double InverseMajorIntegral(const double I) const;

  };// end of class Exponential_Profile
  /*!
    \class Exponential_Profile
    \brief Represents an exponential matter distribution of the hadron 
    in impact parameter space
  */

  class Gaussian_Profile: public Profile_Function_Base {
  private:

    double m_radius;

  public:
  
    // constructor
    Gaussian_Profile(const double radius);

    // member functions
    virtual double Value(const double b) const;
    virtual double MajorValue(const double b) const;

    virtual double MajorIntegral(const double b) const;
    virtual double InverseMajorIntegral(const double I) const;

  };// end of class Gaussian_Profile
  /*!
    \class Gaussian_Profile
    \brief Represents a gaussian matter distribution of the hadron 
    in impact parameter space
  */

  class Double_Gaussian_Profile: public Profile_Function_Base {
  private:

    double m_radius[2], m_partition, m_rmin, m_rmax;

  public:
  
    // constructor
    Double_Gaussian_Profile(const double radius1,const double radius2,
			    const double partition);

    // member functions
    virtual double Value(const double b) const;
    virtual double MajorValue(const double b) const;

    virtual double MajorIntegral(const double b) const;
    virtual double InverseMajorIntegral(const double I) const;

  };// end of class Double_Gaussian_Profile
  /*!
    \class Double_Gaussian_Profile
    \brief Represents a double gaussian matter distribution of the hadron 
    in impact parameter space. The core size and the matter fraction of the 
    core are adjustable.
  */

}// end of namespace AMISIC

#endif
