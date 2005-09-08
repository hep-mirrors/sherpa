#ifndef NLL_Single_Sudakov_H
#define NLL_Single_Sudakov_H

#include "NLL_Branching_Probability_Base.H" 
#include "NLL_Sudakov_Base.H"
#include "Gauss_Integrator.H"

namespace SHERPA {

  class NLL_Single_Sudakov : public NLL_Sudakov_Base  {
  private:

    Sudakov::code m_calcmode, m_cutmode;

    NLL_Branching_Probability_Base *p_bp;

    ATOOLS::Gauss_Integrator m_gauss;

  public:

    NLL_Single_Sudakov(NLL_Branching_Probability_Base *,int);

    ~NLL_Single_Sudakov();

    bool   Initialize(double =0.,double =0.);
    double Log(double,double );
    double operator()(double, double);
    double operator()(double);
    double IntGamma(double,double);

  };

  /*!
    \file 
    \brief contains the class SHERPA::NLL_Single_Sudakov
  */
  /*!
    \class NLL_Single_Sudakov
    \brief provides the sudakov form factor for a single given 
           integrated splitting function
  */

}

#endif
