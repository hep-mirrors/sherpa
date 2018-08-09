#ifndef CFPSHOWER_Shower_Kernels_H
#define CFPSHOWER_Shower_Kernels_H

#include "CFPSHOWER++/Shower/Kernel.H"
#include <vector>

namespace CFPSHOWER {
  class Kernels : public CFPSHOWER::Kernel_Vector {
  private:
    double              m_sum;
    std::vector<double> m_integrals;
    Kernel            * p_selected;
  public:
    Kernels() : p_selected(NULL) {};
    ~Kernels();

    double CalcIntegrals(Splitting & split,const ATOOLS::Mass_Selector * ms);
    bool   SelectOne();

    inline const double & Integral()    const { return m_sum; }
    inline Kernel       * GetSelected() const { return p_selected; }
  };

  std::ostream & operator<<(std::ostream &s,Kernels & kernels);
}

#endif
    
