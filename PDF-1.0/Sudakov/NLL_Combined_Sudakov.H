#ifndef NLL_Combined_Sudakov_H
#define NLL_Combined_Sudakov_H

#include "NLL_Sudakov_Base.H"
#include "NLL_Single_Sudakov.H"
#include <vector>

namespace SHERPA {
  class NLL_Combined_Sudakov : public NLL_Sudakov_Base {
    double m_qmax, m_qmin;
    Sudakov::code m_calcmode;
    Sudakov::code m_cutmode;
    double m_qlimit;
    std::vector<NLL_Single_Sudakov*> m_suds;
  public:
    NLL_Combined_Sudakov(int);
    void Add(NLL_Single_Sudakov*);
    void AssignKey(ATOOLS::Integration_Info *const info) {};

    double operator()(double q, double Q);
    double operator()(double);
    double Log(double,double);

    double Qmax();
    double Qmin();
    void SetQmax(double);
    void SetQmin(double);

    ~NLL_Combined_Sudakov();
  };

  inline double NLL_Combined_Sudakov::Qmax() { return m_qmax; }
  inline double NLL_Combined_Sudakov::Qmin() { return m_qmin; }
  inline void NLL_Combined_Sudakov::SetQmax(double q) {
    m_qmax = q;
  }
  inline void NLL_Combined_Sudakov::SetQmin(double q) {
    m_qmin = q;
  }

  /*!
    \file 
    \brief constains the class SHERPA::NLL_Combined_Sudakov
  */
  /*!
    \class NLL_Combined_Sudakov
    \brief provides the sudakov form factor for a sum of integrated splitting functions
  */

}
#endif
