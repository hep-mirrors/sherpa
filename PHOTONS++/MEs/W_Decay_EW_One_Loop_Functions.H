#ifndef PHOTONS_MEs_W_Decay_EW_One_Loop_Functions_H
#define PHOTONS_MEs_W_Decay_EW_One_Loop_Functions_H

#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"

namespace PHOTONS {
  class W_Decay_EW_One_Loop_Functions: public EW_One_Loop_Functions_Base {
  private:
    // Quantities ending in "p" denote value for isospin partner

    // Vector and axial couplings, charges, combinations
    double m_Qf, m_Qfp, m_If, m_Ifp;
    // Masses squared
    double m_m2;
    double m_m2p;
    
    int m_ew;


  public:
    W_Decay_EW_One_Loop_Functions(const double& s, const double& m2, const double& m2p,
				  const double& Qf, const double& Qfp, 
				  const double& If, const double& Ifp, 
				  const double& mu2, const int& ew);
    virtual ~W_Decay_EW_One_Loop_Functions();

    METOOLS::DivArrC FQED(const int& ind, const int& hel);
    METOOLS::DivArrC FAn();
    METOOLS::DivArrC FZn();
    METOOLS::DivArrC FZa();
    METOOLS::DivArrC CT_L();
    METOOLS::DivArrC CT_R();
    METOOLS::DivArrC V_A(const Complex& m2);
    METOOLS::DivArrC V_B_minus(const Complex& m12, const Complex& m22);
  };
}
    

#endif
