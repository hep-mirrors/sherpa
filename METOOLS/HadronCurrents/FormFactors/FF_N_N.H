#ifndef METOOLS_Current_Library_FF_N_N_H
#define METOOLS_Current_Library_FF_N_N_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// Nucleon form factors for the quasi-elastic N -> N' current VA_N_N,
  /// Tab. "Form factors used in the vertex rules" of implementation_1.tex.
  ///
  /// FORM-FACTOR-indexed, like FF_P_X: one object per invariant function,
  /// selected by its name,
  ///   F1p, F2p, F1n, F2n : Dirac/Pauli form factors of p and n,
  ///   FA                 : nucleon axial form factor F_A.
  /// The boson-dependent combinations (F^V = F^p - F^n, the Z mixtures)
  /// and the PCAC F_P are assembled in VA_N_N, which knows the boson.
  ///
  /// The model number (FORM_FACTOR slot of FF_Parameters) is
  ///   0  : zero
  ///   1  : pointlike - F1p=1, F1n=F2p=F2n=0, F_A=g_A
  ///   10 : Kelly Sachs parametrisation           (F1p..F2n, default)
  ///   20 : Meyer et al. deuterium z-expansion    (FA, default)
  ///   21 : lattice-QCD average z-expansion       (FA)
  ///
  /// Q^2 = -(p_0 - p_1)^2 with m_pi[0] the outgoing, m_pi[1] the incoming
  /// nucleon. m_N and m_pi are Sherpa's own flavour masses.
  ///
  /// Parameters (defaults from the tex input table):
  ///   MU_P 2.793, MU_N -1.913, LAMBDA2 0.7174,
  ///   KELLY_AN 1.7, KELLY_BN 3.3, GA 1.2754.
  //////////////////////////////////////////////////////////////////////////
  class FF_N_N : public FormFactor_Base {
  public:
    enum class Model : int {
      None = 0, Point = 1,
      Kelly = 10,
      Meyer = 20,
      Lattice = 21
    };
  private:
    enum class Which { F1p, F2p, F1n, F2n, FA };
    Which  m_which;
    double m_mN, m_muP, m_muN, m_L2, m_An, m_Bn, m_gA;
    double m_tc, m_t0, m_FA0;
    std::vector<double> m_ak;

    double Kelly(const double & a1,const double & b1,const double & b2,
                 const double & b3,const double & tau) const;
    double Z(const double & Q2) const;
    double ZSeries(const double & Q2) const;
  public:
    FF_N_N(const FF_Parameters & params);
    ~FF_N_N() {}
    double  Value(const double & Q2) const;
    Complex operator()(const ATOOLS::Vec4D_Vector& moms) override;
  };
}

#endif
