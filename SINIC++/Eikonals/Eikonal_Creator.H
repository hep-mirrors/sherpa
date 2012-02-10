#ifndef SINIC_Eikonals_Eikonal_Creator_H
#define SINIC_Eikonals_Eikonal_Creator_H

#include "SINIC++/Eikonals/Omega_ik.H"
#include "SINIC++/Eikonals/Form_Factors.H"

namespace SINIC {
  class Analytic_Contributor : public ATOOLS::Function_Base {
  private:
    Form_Factor * p_ff;
    double        m_Delta, m_Y;
    int           m_sign;
  public:
    Analytic_Contributor(Form_Factor * ff,const double & Delta,const double & Y,
			 const int & sign) :
      p_ff(ff), m_Delta(Delta), m_Y(Y), m_sign(sign) {}
    double operator()(const double & b,const double & y) const;
  };

  class Analytic_Eikonal : public ATOOLS::Function_Base {
  private:
    double m_Delta, m_Y, m_kappa_i, m_kappa_k, m_beta02, m_Lambda2, m_prefactor;
  public:
    Analytic_Eikonal(const double & Delta,const double & Y,
		     const double & kappa_i,const double & kappa_k,
		     const double & beta02,const double & Lambda2) :
      m_Delta(Delta), m_Y(Y), m_kappa_i(kappa_i), m_kappa_k(kappa_k),
      m_beta02(beta02), m_Lambda2(Lambda2), 
      m_prefactor(m_beta02*m_Lambda2/(4.*M_PI)*exp(2.*m_Delta*m_Y)*
		  (1.+m_kappa_i)*(1.+m_kappa_k)/(2.+m_kappa_i+m_kappa_k)) {}
    double operator()(const double & B) const;
  };

  /*!
    \class Eikonal_Creator
    \brief Initialises and fills the single channel eikonals, \f$\Omega_{ik}\f$,
    represented by SINIC::Omega_ik.

    To this end, the eikonal creator solves the coupled differential 
    equations for the two terms \f$Omega_{i(k)}\f$ and \f$\Omega_{(i)k}\f$, individually 
    implemented in the class SINIC::Eikonal_Contributor.  Having them at hand, this class
    also convolutes them by integrating out the individual impact parameters
    \f$\vec b_\perp^{(1,2)}\f$.
  */
  class Eikonal_Creator {
  private:
    Form_Factor     * p_ff1, * p_ff2;
    double            m_lambda, m_Delta, m_beta02;
    absorption::code  m_absorp;
    double            m_originalY, m_cutoffY, m_Y;
    double            m_Bmin, m_Bmax;
    int               m_Bsteps;    
    int               m_test;
    double            m_accu;
    double            m_b1min, m_b2min, m_b1max, m_b2max;
    double            m_ff1max, m_ff2max; 
    int               m_ff1steps, m_ff2steps;

    void   CreateEikonalTerms(Omega_ik * Omega_ik);  
    void   CreateImpactParameterGrid(Omega_ik * Omega_ik);
    int    FixBorders(DEQ_Solver * solver,
		      const double & ff1,const double & ff2,
		      const int & steps=0);

    void   TestEikonal(Omega_ik * omegaik) const;
  public:
    Eikonal_Creator() {}
    Eikonal_Creator(const int & m_test);
    inline ~Eikonal_Creator() {}    

    Omega_ik * CreateEikonal(Form_Factor * ff1,Form_Factor * ff2);
  };
}

#endif
