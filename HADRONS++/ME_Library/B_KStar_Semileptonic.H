#ifndef HADRONS_ME_Library_B_KStar_Semileptonic_H
#define HADRONS_ME_Library_B_KStar_Semileptonic_H

#include "HADRONS++/ME_Library/HD_ME_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Math/MyComplex.H"

namespace HADRONS {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class B_KStar_Semileptonic : public HD_ME_Base {
    // indices: 0=B, 1=K*, 2=lepton-, 3=lepton+
    Complex m_global;
    bool   m_LD;
    int    m_iB, m_iK, m_ilm, m_ilp;
    double m_mB, m_mKhat, m_mbhat, m_mc, m_ms, m_mb;
    double m_C1,m_C2,m_C3,m_C4,m_C5,m_C6,m_C7eff,m_C9,m_C10;
    double m_A1_0,m_A2_0,m_A0_0,m_V_0,m_T1_0,m_T2_0,m_T3_0;
    double m_A1_c1,m_A2_c1,m_A0_c1,m_V_c1,m_T1_c1,m_T2_c1,m_T3_c1;
    double m_A1_c2,m_A2_c2,m_A0_c2,m_V_c2,m_T1_c2,m_T2_c2,m_T3_c2;
    double m_A1_c3,m_A2_c3,m_A0_c3,m_V_c3,m_T1_c3,m_T2_c3,m_T3_c3;
    double m_Vub, m_Vtb;

    public:
    B_KStar_Semileptonic(ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name):
      HD_ME_Base(flavs,n,decayindices,name) {};
    void   operator()(
	const ATOOLS::Vec4D  * _p,
	METOOLS::Spin_Amplitudes     * amps);
    void   SetModelParameters( GeneralModel _md );
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
