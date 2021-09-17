#ifndef CFPSHOWER_Shower_Kernels_H
#define CFPSHOWER_Shower_Kernels_H

#include "CFPSHOWER++/Shower/Kernel.H"
#include <vector>

namespace CFPSHOWER {
  class Kernels : public CFPSHOWER::Kernel_Vector {
  private:
    std::vector<double> m_integrals;
  public:
    Kernels() {}
    ~Kernels();

    double CalcIntegrals(Splitting & split,Configuration & config,
			 const ATOOLS::Mass_Selector * ms);
    inline std::vector<double> & GetIntegrals() { return m_integrals; }
  };

  std::ostream & operator<<(std::ostream &s,Kernels & kernels);
}

#endif
    
