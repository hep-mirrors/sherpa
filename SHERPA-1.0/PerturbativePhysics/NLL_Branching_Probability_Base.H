#ifndef NLL_Branching_Probability_Base_H
#define NLL_Branching_Probability_Base_H

#include "Function_Base.H"

namespace SHERPA {
  class NLL_Branching_Probability_Base : public ATOOLS::Function_Base {
    double m_qmax, m_qmin;
  public:
    virtual double Gamma(double q, double Q) = 0;
    virtual double IntGamma(double q, double Q) = 0;
    virtual double operator()(double q, double Q);
    virtual double operator()(double);

    double Qmax();
    double Qmin();
    void SetQmax(double);
    void SetQmin(double);
  };

  inline double NLL_Branching_Probability_Base::Qmax() { return m_qmax; }
  inline double NLL_Branching_Probability_Base::Qmin() { return m_qmin; }
  inline void NLL_Branching_Probability_Base::SetQmax(double q) {
    m_qmax = q;
  }
  inline void NLL_Branching_Probability_Base::SetQmin(double q) {
    m_qmin = q;
  }

}


#endif
