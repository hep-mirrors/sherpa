#ifndef NLL_Branching_Probability_Base_H
#define NLL_Branching_Probability_Base_H

#include "Function_Base.H"

namespace SHERPA {

  class NLL_Branching_Probability_Base: public ATOOLS::Function_Base {
  protected:

    std::string m_name;
    double m_qmin,m_qmax;

  public:

    virtual double Gamma(double q, double Q) = NULL;
    virtual double IntGamma(double q, double Q) = NULL;

    virtual double operator()(double q, double Q=1000.);
    virtual double operator()(double q);
    virtual double operator()();

    inline void SetQmin(double q) { m_qmin = q; }
    inline void SetQmax(double q) { m_qmax = q; }

    inline std::string Name() const { return m_name; }

    inline double Qmin() const { return m_qmin; }
    inline double Qmax() const { return m_qmax; }

  };

  /*!
    \file 
    \brief Contains the class SHERPA::NLL_Branching_Probability_Base
  */
  /*!
    \class NLL_Branching_Probability_Base
    \brief This class represents a prototype for a branching probability 
           (integrated splitting function), which can be used in the evaluation of
           Sudakov form factors (cf. class NLL_Sudakov)
    
    A single integrated splitting functions has the form
    \f[ 
         \Gamma(Q,q) = \int_{z-(q/Q}^{z+(q/Q)} dz \frac{\alpha_S(q)}{\Pi} P(z)\,,
    \f]
    where \f$\alpha_S\f$ is the (running) strong coupling and  \f$P(z)\f$ is the
    splitting kernel.

    The class provides methods to access
     - the braching probability \f$\Gamma(Q,q)\f$ via Gamma()
     - as well as its integral via IntGamma(), which is used
       as the basis of Sudakov form factors.
    .
  */
  /*!
    \fn  virtual double NLL_Branching_Probability_Base::Gamma(double q, double Q)
    \brief returns the value of a specific braching probability.

    Returns the value of a specific braching probability:
    \f[ 
         \Gamma(Q,q) = \int_{z^-(q/Q}^{z^+(q/Q)} dz \frac{\alpha_S(q)}{\Pi} P(z)\,.
    \f]    
  */
  /*!
    \fn double NLL_Branching_Probability_Base::operator()(double q) 
    \brief returns the value of specific braching probability with a predefined
     upper scale Q (cf. SetQmax()).
  */
  /*!
    \fn double NLL_Branching_Probability_Base::operator()(double q, double Q) 
    \brief alternative access method to Gamma()
  */
  /*!
   \fn    virtual double NLL_Branching_Probability_Base::IntGamma(double Q0, double Q)
   \brief provides the integrated braching probability

   Returns the value of a integrated braching probability:
   \f[
      \int_{Q_0}^{Q} \,dq \, \Gamma(Q,q)\,.
   \f]
  */
  /*!
    \fn     double NLL_Branching_Probability_Base::Qmax();
    \brief  returns the fixed upper scale
  */
  /*!
    \fn    double NLL_Branching_Probability_Base::Qmin();
    \brief returns the fixed lower scale
  */
  /*! 
    \fn void NLL_Branching_Probability_Base::SetQmax(double);
    \brief defines the upper scale used by the operator()(double)
  */
  /*!
    \fn void NLL_Branching_Probability_Base::SetQmin(double);
    \brief defines the fixed lower scale 
  */

}

#endif
