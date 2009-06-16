#ifndef AMISIC_Model_Reggeon_Trajectory_H
#define AMISIC_Model_Reggeon_Trajectory_H
/*!
  \file Reggeon_Trajectory.H
  \brief Declares the class Reggeon_Trajectory
*/

namespace AMISIC {

  class Reggeon_Trajectory {
  private:
    
    double m_alpha_0, m_alpha_prime, m_s, m_s_0;
    
  public:
    
    // constructor
    Reggeon_Trajectory(const double alpha_0,const double alpha_prime);
    /*!
      \fn Reggeon_Trajectory(const double alpha_0,const double alpha_prime)
      \brief The standard constructor. Initializes the Reggeon_Trajectory
      with a given offset and a given slope.
    */
    
    // destructor
    ~Reggeon_Trajectory();
    
    // member functions
    void Fit(const double t,const double dsigma);
    /*!
      \fn void Fit(const double t,const double dsigma)
      \brief Fits the free parameter, \f$s_0\f$ of the Reggeon_Trajectory 
      to the given cross section. 
    */
    
    double DSigma(const double t) const;
    /*!
      \fn double double DSigma(const double t) const
      \brief Returns the differential cross section.
    */
    double GetT(const double tmin,const double tmax,const double ran) const;
    /*!
      \fn double GetT(const double tmin,const double tmax,
      const double ran) const
      \brief Generates an invariant mass squared.
    */
    
    // inline functions
    inline double operator()(const double t=0.0) const { return DSigma(t); }

    inline void SetS(const double s) { m_s=s; }
    
  };// end of class Reggeon_Trajectory
  /*!
    \class Reggeon_Trajectory
    \brief Represents a Reggeon trajectory
    
    This class is employed to generate invariant masses squared, \f$t\f$
    according to reggeon trajectories. Therefore the common form of the
    trajectory, i.e.
    \f[
    \frac{{\rm d}\sigma}{{\rm d}t}\,=\;\left(\frac{s}{s_0}\right)
    ^{2\alpha_0-2\alpha't-2}
    \f]
    is used.
  */

}// end of namespace AMISIC

#endif
