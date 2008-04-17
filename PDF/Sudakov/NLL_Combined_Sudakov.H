#ifndef NLL_Combined_Sudakov_H
#define NLL_Combined_Sudakov_H

#include "NLL_Sudakov_Base.H"
#include "NLL_Single_Sudakov.H"
#include <vector>

namespace SHERPA {

  class NLL_Combined_Sudakov : public NLL_Sudakov_Base {
  private:

    Sudakov::code m_calcmode;
    Sudakov::code m_cutmode;
    std::vector<NLL_Single_Sudakov*> m_suds;

  public:

    NLL_Combined_Sudakov(int);

    ~NLL_Combined_Sudakov();

    void Add(NLL_Single_Sudakov*);

    double operator()(double q, double Q);
    double operator()(double);
    double Log(double,double);

  };
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
