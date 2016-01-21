#ifndef SHRIMPS_Eikonals_Eikonal_Creator_H
#define SHRIMPS_Eikonals_Eikonal_Creator_H

#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Eikonals/Form_Factors.H"

namespace SHRIMPS {
  /*!
    \class Eikonal_Creator
    \brief Initialises and fills the single channel eikonals, 
    \f$\Omega_{ik}\f$, represented by SHRIMPS::Omega_ik.

    To this end, the eikonal creator solves the coupled differential 
    equations for the two terms \f$Omega_{i(k)}\f$ and \f$\Omega_{(i)k}\f$, 
    individually implemented in the class SHRIMPS::Eikonal_Contributor.  
    Having them at hand, this class also convolutes them by integrating out 
    the individual impact parameters \f$\vec b_\perp^{(1,2)}\f$.
  */

  class Eikonal_Creator {
  private:
    Form_Factor      * p_ff1, * p_ff2;
    Eikonal_Parameters m_params;
    int                m_Bsteps,m_ff1steps, m_ff2steps; 

    void FillBYGrids(Eikonal_Contributor * omegai,
		       Eikonal_Contributor * omegak);
    void CreateImpactParameterGrid(Omega_ik * Omega_ik);
    void FixGridAndBorders(DEQ_Solver * solver,int & steps,
			   const double & ff1,const double & ff2);

  public:
    Eikonal_Creator();
    ~Eikonal_Creator() {}    

    void       SetFormFactors(Form_Factor * ff1,Form_Factor * ff2);
    Omega_ik * InitialiseEikonal();
  };
}

#endif
