#ifndef SINIC_Eikonals_Omega_ik_H
#define SINIC_Eikonals_Omega_ik_H

#include "SINIC++/Eikonals/Form_Factors.H"
#include "SINIC++/Eikonals/Eikonal_Contributor.H"
#include "SINIC++/Tools/DEQ_Solver.H"
#include "SINIC++/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"

namespace SINIC {
  /*!
    \class Omega_ik
    \brief The single channel eikonal \f$\Omega_{ik}\f$ contributing to the
    various cross sections.  

    The single channel eikonal \f$\Omega_{ik}\f$ contributing to the
    various cross sections.  It is constructed from the individual terms
    \f$\Omega_{i(k)}\f$ and \f$\Omega_{k(i)}\f$ by integration over the
    relative impact parameters.  This is something that needs to be done,
    my feeling is that we should leave this to some external class
    SINIC::Eikonal::Creator.

    Therefore, this class should know its two SINIC::Eikonal_Contributors, 
    which represent the two terms \f$\Omega_{i(k)}\f$ and \f$\Omega_{k(i)}\f$.  
    It should also have access to the two form factors \f$F_{i,k}\f$ for the 
    actual simulation. For the cross sections, it needs an operator yielding
    \f$\Omega_{ik}(Y,\,y,\,B_\perp)\f$, I could imagine, however, that the 
    depedence on the parameter \f$Y\f$ can be made obsolete.  For checking 
    purposes, though, I guess a dependence on the rapidity point \f$y\f$ 
    would be good to maintain: in principle \f$\Omega_{ik}(Y,\,y,\,B_\perp)\f$ 
    should not depend on it.
  */
  class Omega_ik : public ATOOLS::Function_Base {
  private:
    Form_Factor         * p_ff1, * p_ff2;
    double                m_lambda, m_Delta;
    absorption::code      m_absorp;
    double                m_originalY, m_cutoffY, m_Y, m_Ysteps;
    Eikonal_Contributor   m_Omegaik, m_Omegaki;
    double                m_Bmin, m_Bmax, m_deltaB;
    double                m_singletwt;
    double                m_sigmaInelastic;
    std::vector<double>   m_gridB, m_gridBmax, m_gridD;
    double                m_eta, m_prefactor;
    int                   m_test;
    double                m_diffmax;

  public:
    Omega_ik(Form_Factor * ff1,Form_Factor * ff2,
	     const int & Bsteps,const int & m_test=0);
    ~Omega_ik();

    Eikonal_Contributor * GetSingleTerm(const int & i=0);

    std::vector<double> * GetImpactParameterGrid()        { return &m_gridB; }
    std::vector<double> * GetImpactParameterMaximumGrid() { return &m_gridBmax; }
    void SetImpactParameterInterval(const double & Bmin,const double & Bmax) {
      m_Bmin   = Bmin; 
      m_Bmax   = Bmax; 
      m_deltaB = (m_Bmax-m_Bmin)/double(m_gridB.size()-1);
    }

    double operator()(const double & B) const;
    double Maximum(const double & B) const;

    ATOOLS::Vec4D SelectB1B2(double & b1,double & b2,const double & B);
    double MaximalEmissionProbability(const double & b1,const double & b2);
    double EmissionWeight(const double & b1,const double & b2,const double & y,
			  const double & sup=1.);
    double SingletWeight(const double & b1,const double & b2,
			 const double & y1,const double & y2,
			 const double & sup=1.,const int & nbeam=1);
    double OctetWeight(const double & b1,const double & b2,
		       const double & y1,const double & y2,
		       const double & sup=1.,const int & nbeam=1);
    double Sum(const double & b1,const double & b2,const double & y);
    double RescatterProbability(const double & b1,const double & b2,
				const double & y1,const double & y2,
				const double & sup=1.,const int & nbeam=0);

    void   PrepareQT(const double & b1,const double & b2);

    const double   Prefactor() const { 
      return ATOOLS::sqr(p_ff1->Prefactor()*p_ff2->Prefactor()); 
    }
    const double & lambda() const    { return m_lambda; }
    const double & Delta() const     { return m_Delta; }
    const double & Kappa_i() const   { return p_ff1->Kappa(); }
    const double & Kappa_k() const   { return p_ff2->Kappa(); }
    const double & Lambda2() const   { return p_ff1->Lambda2(); }
    const double & Xi() const        { return p_ff1->Xi(); }
    const double & Ymax() const      { return m_Y; }
    double EffectiveIntercept(double b1=-1.,double b2=-1.,const double & y=0.);

    Form_Factor * FF1() const  { return p_ff1; }
    Form_Factor * FF2() const  { return p_ff2; }

    void  SetSigmaInelastic(const double & sigma) { m_sigmaInelastic = sigma; }
    const double & Sigma_Inelastic() const        { return m_sigmaInelastic; }
  };
}

#endif
