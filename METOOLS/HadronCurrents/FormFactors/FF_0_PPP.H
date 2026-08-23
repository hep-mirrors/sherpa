#ifndef METOOLS_Current_Library_FF_0_PPP_H
#define METOOLS_Current_Library_FF_0_PPP_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  enum FF_0_PPP_mode {
    piP_pi0_pi0       = 1,
    piM_piP_piP       = 2,
    // K^+ K^- pi^- (cf. VA_0_KKPi.H for the fixed 0=K+,1=K-,2=pi-
    // ordering convention this mode assumes).
    KP_KM_piM         = 3,
    // pi^0 pi^0 K^- (Bose-symmetric like piP_pi0_pi0/piM_piP_piP above -
    // reuses VA_0_PiPiPi's Current class and its p1=q1,p2=q2,p3=q3
    // convention). Finkemeier-Mirkes hep-ph/9503474 (FM95), Tabs.I/II.
    pi0_pi0_KM        = 4,
    // Below: new "KS"-mode (Kuehn-Santamaria-style) channels from
    // Finkemeier & Mirkes, hep-ph/9503474 (FM95) and hep-ph/9601275
    // (FM96, scalar Kpi only - see FF_0_PP.C for that one). None of
    // these are Bose-symmetric (all three mesons distinct, or KS/KL
    // treated as their own distinct flavours), so they are wired to a
    // new pair of generic Current classes (VA_0_KPiK, VA_0_KPiPi; see
    // those .H/.C files) rather than reusing VA_0_PiPiPi/VA_0_KKPi.
    // FixMode assigns p1=q1,p2=q2,p3=q3 in FM95's own labeling for
    // each channel (documented per-channel at the FixMode() site and
    // repeated at each FF class), i.e. the SAME base-class convention
    // as piP_pi0_pi0/piM_piP_piP above: s1_arg=(p1+p2)^2=paper's s3,
    // s2_arg=(p1+p3)^2=paper's s2, and paper's s1=(p2+p3)^2 must be
    // reconstructed via momentum conservation.
    KM_piM_KP         = 5,  // K^-(q1) pi^-(q2) K^+(q3)   [FM95 Tab.I/II row 1]
    K0_piM_K0bar      = 6,  // K^0(q1) pi^-(q2) K0bar(q3) [row 2]
    KS_piM_KS         = 7,  // K_S(q1) pi^-(q2) K_S(q3)   [row 3]
    KS_piM_KL         = 8,  // K_S(q1) pi^-(q2) K_L(q3)   [row 4]
    KL_piM_KL         = 9,  // K_L(q1) pi^-(q2) K_L(q3)   [Eq.(50), = -KSpiKS-like G's]
    KM_pi0_K0         = 10, // K^-(q1) pi^0(q2) K^0(q3)   [row 5]
    KM_piM_piP        = 11, // K^-(q1) pi^-(q2) pi^+(q3)  [row 7]
    piM_K0bar_pi0     = 12, // pi^-(q1) K0bar(q2) pi^0(q3)[row 8]
    // eta(q1) pi^-(q2) pi^0(q3) / eta'(q1) pi^-(q2) pi^0(q3):
    // tau_two_meson_currents_KS_RChiT.tex Sec."eta pi pi and eta' pi
    // pi": G-parity kills the axial form factors identically in the
    // isospin limit - only the anomalous vector (F5/FS) current
    // survives (Eq.(eq:etapipi-current)). Same p1=q1,p2=q2,p3=q3
    // convention as above (q1=eta(or eta'), q2=pi^-, q3=pi^0, per the
    // note's own labeling) - NOT interchangeable (unlike the Keta/Kpi
    // two-meson case) since eta(') is not Bose-symmetric with either
    // pion, so FixMode requires this exact order.
    EtaPiPi_pi0       = 13,
    EtaprimePiPi_pi0  = 14,
    unknown           = 999
  };
  
  class FF_0_PPP_Base : public FormFactor_Base {
  protected:
    FF_0_PPP_mode       m_mode;
    Complex             m_norm;
    // True if the ACTUAL final-state kaon in piM_K0bar_pi0 is K_S or
    // K_L rather than the bare (unobservable) K^0/K0bar flavour
    // eigenstate - same reasoning as FF_0_PP.H's m_isKSKL (K_S/K_L are
    // strangeness-eigenstate MIXTURES, picking up an ADDITIONAL
    // 1/sqrt(2) amplitude factor on top of the isospin Clebsch already
    // in m_norm for this channel). Not relevant for the KpiK family's
    // KS_piM_KS/KS_piM_KL/KL_piM_KL modes, which already recognise
    // K_S/K_L directly with their own dedicated normalisation (Eq.50).
    bool                m_isKSKL;

    void    FixMode();
    virtual Complex FF_KS(const double & s123,const double & s1,const double & s2);
    virtual Complex FF_RChiPT(const double & s123,const double & s1,const double & s2);
    // Base-class arg convention (matches FF_KS/FF_RChiPT above, fixed by
    // FF_0_PPP_Base::operator()): s123=q^2=(p1+p2+p3)^2, s1=(p1+p2)^2,
    // s2=(p1+p3)^2. NOTE this differs from the s1,s2 NAMING used in
    // 1203.3955/1509.09140 (there, s1=(p2+p3)^2, s2=(p1+p3)^2,
    // s3=(p1+p2)^2) - i.e. this class's "s1" argument is the PAPERS'
    // s3, and the papers' s1 must be reconstructed as
    // q^2 - s2_arg - s1_arg + m1^2+m2^2+m3^2. Each RChL2012 override
    // does this conversion explicitly at the top of the function to
    // avoid silent mislabeling.
    virtual Complex FF_RChL2012(const double & s123,const double & s1,const double & s2) {
      return Complex(0.,0.);
    }
  public:
    FF_0_PPP_Base(const FF_Parameters & params); 
    ~FF_0_PPP_Base() {}

    Complex operator()(const ATOOLS::Vec4D_Vector& moms);
    inline const Complex & GetNorm() const { return m_norm; }
  };
}

#endif
