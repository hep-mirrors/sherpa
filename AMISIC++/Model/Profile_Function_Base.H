#ifndef AMISIC_Model_Profile_Function_Base_H
#define AMISIC_Model_Profile_Function_Base_H
/*!
  \file Profile_Function_Base.H
  \brief Declares the class Profile_Function_Base
*/

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Org/Message.H"
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
  /*!
    \struct pft
    \brief Declares abbreviations for profile function names.

    This struct is employed to declare abbreviations 
    for profile function names.
  */

  std::ostream &operator<<(std::ostream &ostr,const pft::code code);

  class Profile_Function_Base;

  class Differential_Overlap: public ATOOLS::Function_Base {
  private:

    Profile_Function_Base *p_owner;

  public:
    
    // constructor
    Differential_Overlap(Profile_Function_Base *const owner);
    /*!
      \fn Differential_Overlap(Profile_Function_Base *const owner)
      \brief The standard constructor.
    */
    
    // member functions
    double operator()(const double b);
    /*!
      \fn double operator()(const double b)
      \brief Returns \f$\tilde{O}(b)\f$.
    */
    
  };// end of class Differential_Overlap
  /*!
    \class Differential_Overlap
    \brief Employed to determine the differential overlap between 
    two colliding hadrons in impact parameter space 

    The class Differential_Overlap is employed to calculate the 
    differential overlap of two hadron profiles in impact parameter 
    space, which is given by
    \f[ 
    \tilde{O}(b)\,=\;\int{\rm d}^3{\bf x}\,{\rm d} t\,
    \rho_{1\rm b}\left(x-{\textstyle\frac{1}{2}}\,b,y,z-vt\right)\,
    \rho_{2\rm b}\left(x+{\textstyle\frac{1}{2}}\,b,y,z+vt\right)\;.
    \f]    
  */

  class Interaction_Probability: public ATOOLS::Function_Base {
  private:

    Profile_Function_Base *p_owner;

  public:
    
    // constructor
    Interaction_Probability(Profile_Function_Base *const owner);
    /*!
      \fn Interaction_Probability(Profile_Function_Base *const owner)
      \brief The standard constructor.
    */
    
    // member functions
    double operator()(const double b);
    /*!
      \fn double operator()(const double b)
      \brief Returns \f$P_{\rm int}(b)\f$.
    */
    
  };// end of class Interaction_Probability
  /*!
    \class Interaction_Probability
    \brief Employed to determine interaction probabilities.

    The class Interaction_Probability is employed to determine the 
    probability for an interaction between two colliding hadrons with
    the differential overlap \f$\tilde{O}(b)\f$ in impact parameter space, 
    which according to Poissonian statistics is given by
    \f[ 
    P_{\rm int}(b)\,=\;1-\exp\left\{\,-\left<\tilde{n}(b)\right>\,\right\}\;.
    \f]    
  */
  
  class Profile_Function_Base {
  protected:
    
    Differential_Overlap    *p_overlap;
    Interaction_Probability *p_probability;

    pft::code m_type;

    double m_bmin, m_bmax, m_omin, m_omax;
    double m_kfactor, m_omean, m_norm;

    // constructor
    Profile_Function_Base(const pft::code code,
			  const double bmin,const double bmax);

    template <class Profile_Type>
    static Profile_Function_Base *
    CreateProfile(const std::string &type,
		  const std::vector<double> &parameters);

  public:
  
    // destructor
    virtual ~Profile_Function_Base();

    static Profile_Function_Base *
    SelectProfile(const std::string &type,
		  const std::vector<double> &parameters);
    /*!
      \fn static Profile_Function_Base *
      SelectProfile(const std::string &type,
      const std::vector<double> &parameters)
      \brief Selects a hadron profile according to the given name 
      and the number of given parameters.
    */

    // member functions
    virtual double Value(const double b) const = 0;
    /*!
      \fn virtual double Value(const double b) const
      \brief Returns the differential overlap \f$\tilde{O}(b)\f$.
    */
    virtual double MajorValue(const double b) const = 0;
    /*!
      \fn virtual double MajorValue(const double b) const
      \brief Returns the differential overlap \f$\tilde{O}(b)\f$
      according to the major.
    */

    virtual double MajorIntegral(const double b) const = 0;
    /*!
      \fn virtual double MajorIntegral(const double b) const
      \brief Returns the integral \f$I(b)\,=\;\int_{b_{min}}^{b}
      {\rm d}b\,\tilde{O}(b)\f$ according to the major.
    */
    virtual double InverseMajorIntegral(const double I) const = 0;
    /*!
      \fn virtual double InverseMajorIntegral(const double I) const
      \brief Returns the inverse integral \f$I^{-1}(b)\f$ 
      according to the major.
    */

    bool CalculateOMean(const double ratio);
    /*!
      \fn bool CalculateOMean(const double ratio)
      \brief Calculates the mean overlap \f$\left<\tilde{O}\right>\f$
    */

    double DiceImpactParameter() const;

    // inline functions
    inline pft::code Type() const { return m_type; }

    inline double BMin() const { return m_bmin; }
    inline double BMax() const { return m_bmax; }

    inline double OMin() const { return m_omin; }
    inline double OMax() const { return m_omax; }

    inline double Norm() const { return m_norm; }

    inline double KFactor() const { return m_kfactor; }
    /*!
      \fn inline double KFactor() const
      \brief Returns the k factor \f$k\f$.
    */
    inline double OMean() const   { return m_omean;   }
    /*!
      \fn inline double OMean() const
      \brief Returns the mean Overlap \f$\left<\tilde{O}\right>\f$.
    */

    inline Differential_Overlap *const Differential() const 
    { return p_overlap; }
    inline Interaction_Probability *const Probability() const  
    { return p_probability; }

    inline double operator()(const double b) const 
    { return Value(b); }
    inline double operator[](const double b) const 
    { return MajorValue(b); }

  };// end of class Profile_Function_Base
  /*!
    \class Profile_Function_Base
    \brief Employed to simulate hadron profiles in impact parameter space

    The class Profile_Function_Base is the base class for all hadron 
    profiles to be used in the simulation of underlying events.
    It employs the classes Differential_Overlap and Interaction_Probability
    to determine the mean overlap of two colliding hadrons, which in the 
    Sherpa framework is given by
    \f[
    \left<\tilde{O}\right>\,=\;\frac{\displaystyle\int{\rm d}^2 {\bf b}\,
    \tilde{O}(b)}
    {\displaystyle\int{\rm d}^2 {\bf b}\,P_{\rm int}(b)}\,=\;\frac{1}{f_c}\,
    \frac{\displaystyle\int{\rm d}^2 {\bf b}\,
    \tilde{O}(b)\,P_{\rm int}(b)}
    {\displaystyle\int{\rm d}^2 {\bf b}\,P_{\rm int}(b)}\;.
    \f]
    The mean multiplicity of underlying events then follows from
    \f[ 
    \left<\tilde{n}(b)\right>\,=\;f(b)\left<k\,\tilde{O}\right>
    \quad\quad\quad {\rm where}\quad\quad\quad
    f(b)\,=\;\frac{\tilde{O}(b)}{\left<\tilde{O}\right>}
    \f]
  */

}// end of namespace AMISIC

#endif
