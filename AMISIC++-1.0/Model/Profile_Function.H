#ifndef Profile_Function_H
#define Profile_Function_H

#include "Profile_Function_Base.H"

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

}// end of namespace AMISIC

#endif
