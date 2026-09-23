#ifndef METOOLS_Current_Library_VA_P_PP_H
#define METOOLS_Current_Library_VA_P_PP_H

#include "METOOLS/HadronCurrents/VA_P_X_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/FF_PP_Wave.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// D -> (P1 P2)_{S-wave} l nu : the genuine four-body S-wave current.
  ///
  /// Naming follows the directory convention <in>_<out>: one
  /// pseudoscalar in, two pseudoscalars out. Indices are
  ///     0 = decaying D,  1 = first pseudoscalar,  2 = second,
  /// e.g. for the channel "-321,211,12,-11" (D+ -> K- pi+ e+ nu):
  ///     J2 Indices: [0, 1, 2],   J1 (VA_0_FF) Indices: [3, 4].
  ///
  /// The Lorentz structure is that of VA_P_S with the scalar's
  /// momentum replaced by Q = p_1 + p_2 and its mass by the GENERATED
  /// sqrt(s), multiplied by a complex S-wave lineshape A_S(s):
  ///
  ///   H^mu = i A_S(s) { f_+^S(q^2,s) [ P^mu - (Delta/q^2) q^mu ]
  ///                    + f_0^S(q^2,s) (Delta/q^2) q^mu },
  ///   Q = p_1+p_2,  q = p_D - Q,  P = p_D + Q,  Delta = M^2 - s.
  ///
  /// This is what the quasi-two-body D -> K*_0(1430) entry cannot do.
  /// PDG 2026 measures B(D+ -> (K-pi+)_S e+ nu) = (2.28+-0.11)e-3 and
  /// B(D0 -> (K0bar pi-)_S e+ nu) = (8.5+-0.6)e-4, and near threshold
  /// that S-wave is kappa/K*_0(700)-dominated, i.e. a broad amplitude
  /// with a slowly turning phase rather than a pole at 1.43 GeV.
  ///
  /// !! DOUBLE COUNTING !! this current covers the WHOLE S-wave. A
  /// decay table using it must NOT also carry a quasi-two-body
  /// D -> K*_0(1430) l nu channel, or the S-wave is counted twice.
  /// The shipped Decaydata.yaml drops those four channels.
  ///
  /// INTERFERENCE is now included: S, P and D waves are summed at the
  /// AMPLITUDE level before the current is inserted, so the four-body
  /// angular structure - including the S-P forward-backward asymmetry
  /// that a sum of separate channels cannot produce - comes out of the
  /// construction rather than being imposed.
  ///
  /// The three Lorentz structures are NOT rederived here. They are the
  /// same VA_P_X_Base kernels the quasi-two-body currents call, with
  /// the polarization argument replaced by a wave projector, and with
  /// Lhat normalised so that the narrow-width limit of the P wave
  /// reproduces VA_P_V exactly.
  //////////////////////////////////////////////////////////////////////////
  class VA_P_PP : public VA_P_X_Base {
  private:
    /// S-wave form factors (D -> S structure).
    FF_P_X * p_fplus, * p_fzero;
    /// P-wave form factors (D -> V structure).
    FF_P_X * p_V, * p_A0, * p_A1, * p_A2;
    /// D-wave form factors (D -> T structure).
    FF_P_X * p_VT, * p_A0T, * p_A1T, * p_A2T;
    /// Lineshapes A_L(s), one per wave.
    FF_PP_Wave * p_S, * p_P, * p_D;
    /// Complex wave coefficients.
    Complex m_cS, m_cP, m_cD;
    bool    m_useS, m_useP, m_useD;
    double  m_phase;
    double  m_m1, m_m2;
    /// Per-wave isospin Clebsch-Gordan factors. Derived from the
    /// flavours for K pi; must be supplied for systems whose waves
    /// carry different isospin (pi pi: S is I=0, P is I=1).
    double  m_isoS, m_isoP, m_isoD;
    /// V->PP and T->PP couplings, derived from each propagator's pole
    /// mass and on-shell width. Overridable with G_P / G_D.
    Complex m_gS, m_gP, m_gD;

    /// Clebsch-Gordan coefficient for the given wave. Sets 'known'
    /// false if the pair is not one it can resolve.
    double IsospinCG(const int & wave,bool & known) const;
    double Breakup(const double & s) const;
    Complex VertexCoupling(Propagator_Base * prop,const int & L) const;

    /// Normalised transverse relative momentum, Lhat.Lhat = -1.
    ATOOLS::Vec4C WaveP(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2,
                        const ATOOLS::Vec4D & Q,const double & s) const;
    /// Spin-2 projector built from Lhat, contracted with p_D and /M.
    ATOOLS::Vec4C WaveD(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2,
                        const ATOOLS::Vec4D & Q,const double & s,
                        const ATOOLS::Vec4D & pD) const;
    FF_PP_Wave * MakeWave(const std::string & name,GeneralModel & model);

  public:
    VA_P_PP(const ATOOLS::Flavour_Vector& flavs,
            const std::vector<int>& indices,
            const std::string& name);
    ~VA_P_PP() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
