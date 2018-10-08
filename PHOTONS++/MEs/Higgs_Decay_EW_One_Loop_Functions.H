#ifndef PHOTONS_MEs_Higgs_Decay_EW_One_Loop_Functions_H
#define PHOTONS_MEs_Higgs_Decay_EW_One_Loop_Functions_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"

namespace PHOTONS {
  class Higgs_Decay_EW_One_Loop_Functions: public EW_One_Loop_Functions_Base {
  private:
    // Quantities ending in "p" denote value for isospin partner

    // Vector and axial couplings, charges, combinations
    double m_Qf, m_Qfp, m_If, m_Ifp;
    // Masses squared
    double m_m2, m_m2p;
    // Mass ratios
    double w, mu2w, z, wfp, wf, wh;
    Complex sigma2, v;
    
    Complex f(const int& i);
    Complex h(const int& i);

  public:
    Higgs_Decay_EW_One_Loop_Functions(const double& s, const double& m2, const double& m2p,
				      const double& Qf, const double& If,
				      const double& mu2, const int& ew);
    virtual ~Higgs_Decay_EW_One_Loop_Functions();
    METOOLS::DivArrC FSQED();
    METOOLS::DivArrC FS();
    METOOLS::DivArrC CT_L();
    METOOLS::DivArrC CT_R();

  };
}
    

#endif
