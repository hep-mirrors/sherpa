#ifndef NLL_Single_Sudakov_H
#define NLL_Single_Sudakov_H

#include "Function_Base.H"
#include "MathTools.H"
#include "Gauss_Integrator.H"
#include "Fast_Function.H"

#include "NLL_Branching_Probability_Base.H" 
#include "NLL_Sudakov_Base.H"

namespace SHERPA {

  class NLL_Single_Sudakov : public NLL_Sudakov_Base  {
    double m_qmax, m_qmin;
    Sudakov::code m_calcmode;
    Sudakov::code m_cutmode;
    double m_qlimit;
    NLL_Branching_Probability_Base * p_bp;
    ATOOLS::Gauss_Integrator m_gauss;
    ATOOLS::Fast_Function    m_log_delta;
    std::string m_outpath;
    double m_lastsud;
  public:
    NLL_Single_Sudakov(NLL_Branching_Probability_Base *,int);
    bool Initialize(double =0.,double =0.);
    void   AssignKey(ATOOLS::Integration_Info *const info) {};

    double operator()(double, double);
    double operator()(double);
    double Log(double,double );
    double IntGamma(double,double);
    double Qmax();
    double Qmin();
    void SetQmax(double);
    void SetQmin(double);
    void SetOutPath(const std::string);

    ATOOLS::Gauss_Integrator & GetIntegrator();

    ~NLL_Single_Sudakov();
  };

  inline double NLL_Single_Sudakov::Qmax() { return m_qmax; }
  inline double NLL_Single_Sudakov::Qmin() { return m_qmin; }
  inline void NLL_Single_Sudakov::SetQmin(double q) {
    m_qmin = q;
    p_bp->SetQmin(q);
  }
  inline void NLL_Single_Sudakov::SetQmax(double q) {
    m_qmax = q;
    p_bp->SetQmax(q);
  }
  inline void NLL_Single_Sudakov::SetOutPath(const std::string _m_outpath) {
    m_outpath=_m_outpath;
  }
  inline double NLL_Single_Sudakov::IntGamma(double q,double Q) {
    return  p_bp->IntGamma(q,Q);
  }
  inline ATOOLS::Gauss_Integrator & NLL_Single_Sudakov::GetIntegrator() {
    return m_gauss;
  }

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
