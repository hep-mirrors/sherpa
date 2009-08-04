#ifndef HADRONS_ME_Library_B_K_Semileptonic_H
#define HADRONS_ME_Library_B_K_Semileptonic_H
 
#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS { 

  class B_K_Semileptonic : public HD_ME_Base {
    // indices: 0=B, 1=K, 2=lepton-, 3=lepton+
    double m_global;
    bool m_LD;
    Complex m_cR_T1, m_cR_T2, m_cL_T1, m_cL_T2;
    double m_C1, m_C2, m_C3, m_C4, m_C5, m_C6, m_C7eff, m_C9, m_C10;
    double m_mc, m_mb, m_ms, m_Vtb, m_Vub;
    public:
    B_K_Semileptonic(ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name):
      HD_ME_Base(flavs,n,decayindices,name) {};
    void   operator()(
      const ATOOLS::Vec4D  * _p,
      METOOLS::Spin_Amplitudes     * amps);
    void SetModelParameters( GeneralModel _md );
    Complex C9sehgal(double sHat);
    Complex sehgalld(double sHat);
    Complex g(double shat);
    Complex g0(double shat);
    Complex gc(double shat);
    double Theta(double x);
    Complex gSD(double mhat, double shat);
  };
}

#endif
