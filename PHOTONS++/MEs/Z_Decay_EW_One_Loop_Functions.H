#ifndef PHOTONS_MEs_Z_Decay_EW_One_Loop_Functions_H
#define PHOTONS_MEs_Z_Decay_EW_One_Loop_Functions_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"

namespace PHOTONS {
  class Z_Decay_EW_One_Loop_Functions: public EW_One_Loop_Functions_Base {
  private:
    // Quantities ending in "p" denote value for isospin partner

    // Vector and axial couplings, charges, combinations
    double m_Qf;
    double m_If;
    Complex beta2;
    Complex kappa;
    // Masses squared
    double m_m2;
    double m_m2p;
    // Mass ratios
    Complex wf;
    Complex zf;
    Complex wfp;
    Complex zfp;
    int m_ew;

  public:
    Z_Decay_EW_One_Loop_Functions(const double& s, const double& m2, const double& m2p,
				  const double& Qf, const double& If,
				  const double& mu2, const int& ew);
    virtual ~Z_Decay_EW_One_Loop_Functions();
    METOOLS::DivArrC FAa();
    METOOLS::DivArrC FA1();
    METOOLS::DivArrC FA3();
    METOOLS::DivArrC FV2();
    METOOLS::DivArrC FZa();
    METOOLS::DivArrC FWa();
    METOOLS::DivArrC FWabar();
    METOOLS::DivArrC FWn();
    METOOLS::DivArrC FWnbar();
    METOOLS::DivArrC CT_L();
    METOOLS::DivArrC CT_R();
  };
}
    

#endif
