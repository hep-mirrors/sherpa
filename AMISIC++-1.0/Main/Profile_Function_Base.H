#ifndef Profile_Function_Base_H
#define Profile_Function_Base_H

#include "Function_Base.H"
#include "Message.H"
#include <vector>

namespace AMISIC {

  struct pft {

    enum code {
      none            = 0,
      flat            = 1,
      exponential     = 2,
      gaussian        = 3,
      double_gaussian = 4
    };

  };// end of struct pft

  std::ostream &operator<<(std::ostream &ostr,const pft::code code);

  class Profile_Function_Base;

  class Differential_Overlap: public ATOOLS::Function_Base {
  private:

    Profile_Function_Base *p_owner;

  public:
    
    // constructor
    Differential_Overlap(Profile_Function_Base *const owner);
    
    // member functions
    double operator()(const double b);
    
  };// end of class Differential_Overlap

  class Interaction_Probability: public ATOOLS::Function_Base {
  private:

    Profile_Function_Base *p_owner;

  public:
    
    // constructor
    Interaction_Probability(Profile_Function_Base *const owner);
    
    // member functions
    double operator()(const double b);
    
  };// end of class Interaction_Probability
  
  class Profile_Function_Base {
  protected:
    
    Differential_Overlap    *p_overlap;
    Interaction_Probability *p_probability;

    pft::code m_type;

    double m_bmin, m_bmax, m_omin, m_omax;
    double m_kfactor, m_omean, m_norm;

    // constructor
    Profile_Function_Base(const pft::code code,const double bmin,const double bmax);

    template <class Profile_Type>
    static Profile_Function_Base *CreateProfile(const std::string &type,
						const std::vector<double> &parameters);

  public:
  
    // destructor
    virtual ~Profile_Function_Base();

    static Profile_Function_Base *SelectProfile(const std::string &type,
						const std::vector<double> &parameters);

    // member functions
    virtual double Value(const double b) const = 0;
    virtual double MajorValue(const double b) const = 0;

    virtual double MajorIntegral(const double b) const = 0;
    virtual double InverseMajorIntegral(const double I) const = 0;

    bool CalculateOMean(const double ratio);

    // inline functions
    inline pft::code Type() const { return m_type; }

    inline double BMin() const { return m_bmin; }
    inline double BMax() const { return m_bmax; }

    inline double OMin() const { return m_omin; }
    inline double OMax() const { return m_omax; }

    inline double Norm() const { return m_norm; }

    inline double operator()(const double b) const { return Value(b);      }
    inline double operator[](const double b) const { return MajorValue(b); }

    inline double operator+(const double b) const { return MajorIntegral(b);        }
    inline double operator-(const double I) const { return InverseMajorIntegral(I); }

    inline double KFactor() const { return m_kfactor; }
    inline double OMean() const   { return m_omean;   }

    Differential_Overlap    *const Differential() const { return p_overlap;     }
    Interaction_Probability *const Probability() const  { return p_probability; }

  };// end of class Profile_Function_Base

}// end of namespace AMISIC

#endif
