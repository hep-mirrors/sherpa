#ifndef SHRIMPS_Eikonals_Omega_ik_H
#define SHRIMPS_Eikonals_Omega_ik_H

#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"
#include "SHRiMPS/Eikonals/Analytic_Contributor.H"
#include "SHRiMPS/Eikonals/Eikonal_Weights.H"
#include "SHRiMPS/Tools/DEQ_Solver.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"

namespace SHRIMPS {
  /*!
    \class Omega_ik
    \brief The single channel eikonal \f$\Omega_{ik}\f$ contributing to the
    various cross sections.  

    The single channel eikonal \f$\Omega_{ik}\f$ contributing to the
    various cross sections.  It is constructed from the individual terms
    \f$\Omega_{i(k)}\f$ and \f$\Omega_{k(i)}\f$ by integration over the
    relative impact parameters.  This is something that needs to be done,
    my feeling is that we should leave this to some external class
    SHRIMPS::Eikonal::Creator.

    Therefore, this class should know its two SHRIMPS::Eikonal_Contributors, 
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
    Eikonal_Weights       m_weights;
    Eikonal_Contributor * p_Omegaik, * p_Omegaki;

    double m_bmax;
    double m_Y;
    int    m_Bsteps, m_deltaB, m_Ysteps;
    double m_prefactor,m_sigmaInelastic;
    std::vector<double>   m_gridB, m_gridBmax, m_gridD;   
  public:
    Omega_ik(const Eikonal_Parameters & params);
    ~Omega_ik();

    void SetContributors(Eikonal_Contributor * Omegaik, 
			 Eikonal_Contributor * Omegaki);
    void SetPrefactor(const double & prefactor) { m_prefactor = prefactor; }

    Eikonal_Contributor * GetSingleTerm(const int & i=0);
    Eikonal_Weights     * GetWeights();

    std::vector<double> * GetImpactParameterGrid()        { return &m_gridB; }
    std::vector<double> * GetImpactParameterMaximumGrid() { return &m_gridBmax; }
    
    double operator()(const double & B) const;
    double Maximum(const double & B) const;    
    ATOOLS::Vec4D SelectB1B2(double & b1,double & b2,const double & B);
    
    void   PrepareQT(const double & b1,const double & b2);
    
    double   Prefactor() const { 
      return m_prefactor;
    }
    // const double & lambda() const  { return m_params.lambda; }
    // const double & Delta() const   { return m_params.Delta; }
    // const double & Kappa_i() const { return p_Omegaik->Kappa(); }
    // const double & Kappa_k() const { return p_Omegaki->Kappa(); }
    // const double & Lambda2() const { return p_Omegaik->Lambda2(); }
    // const double & Ymax() const    { return m_params.Ymax; } 
    double EffectiveIntercept(double b1=-1.,double b2=-1.,const double & y=0.);
    
    Form_Factor * FF1() const  { return p_Omegaik->FF(); }
    Form_Factor * FF2() const  { return p_Omegaki->FF(); }
    
    void  SetSigmaInelastic(const double & sigma) { m_sigmaInelastic = sigma; }
    const double & Sigma_Inelastic() const        { return m_sigmaInelastic; }

    void TestIndividualGrids(Analytic_Contributor * ana12,
			     Analytic_Contributor * ana21,
			     const double & Ymax,const std::string & dirname) const;
  };
  
  struct eikcomp{
    bool operator() (const Omega_ik * eik1, const Omega_ik * eik2) const
    {
      if(eik1->FF1()->Number() < eik2->FF1()->Number()) return true;
      else if(eik1->FF1()->Number() > eik2->FF1()->Number()) return false;
      else {
	if (eik1->FF2()->Number() < eik2->FF2()->Number()) return true;
	else return false;
      }
    }
  };

}

#endif
