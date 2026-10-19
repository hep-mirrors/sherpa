#ifndef METOOLS_Current_Library_VA_N_N_H
#define METOOLS_Current_Library_VA_N_N_H

#include "METOOLS/HadronCurrents/Current_Base.H"

namespace METOOLS {
  //////////////////////////////////////////////////////////////////////////
  /// Quasi-elastic nucleon current, Eqs. (VertexNNA), (VertexNNZ),
  /// (VertexNNW) of implementation_1.tex:
  ///
  ///   J^mu = i ubar(p') [ gamma^mu F1 + i sigma^{mu nu} q_nu/(2m_N) F2
  ///                       - ( gamma^mu F_A + q^mu/(2m_N) F_P ) gamma_5 ] u(p)
  ///
  /// with q = p' - p, Q^2 = -q^2, and the boson-dependent form factors
  ///   photon : F_i = F_i^N,                     F_A = F_P = 0
  ///   Z      : F_i = T3 F_i^V - 2 s_W^2 F_i^N,  F_A -> T3 F_A
  ///   W      : F_i = F_i^V, times V_ud
  /// and F^V = F^p - F^n, F_P = 4 m_N^2/(Q^2+m_pi^2) F_A (PCAC).
  ///
  /// The overall i and the coupling convention match VA_0_FF, so that
  /// Exchange_ME's g_V = e^2, (g/2c_W)^2, (g/2sqrt2)^2 complete the vertex.
  /// The form factors carry the nucleon's charges (F1^p(0)=1, F1^n(0)=0),
  /// so v and a are NOT used. The gamma_5 terms are fed to XYZFunc through
  /// the same (cR,cL) = (V+A, V-A) map as VA_0_FF. XYZFunc's cR multiplies
  /// P_R = (1+gamma_5)/2 in both L and Y (hep-ph/0109036, Eq. 24), so the
  /// gamma^mu gamma_5 and the F_P gamma_5 terms share one convention. That
  /// map gives V + A gamma_5, the parity mirror of the tex vertex, exactly
  /// as Exchange_ME's DefaultVA (a = T3) mirrors the Z lepton vertex; the
  /// relative sign, and so every unpolarised rate, is that of the tex.
  /// Give the lepton line of a W exchange v: 1, a: 1 to stay consistent.
  ///
  /// sigma^{mu nu} is removed by the on-shell Gordon identity
  ///   ubar(p') i sigma^{mu nu} q_nu u(p) = ubar(p') [(m'+m) gamma^mu - P^mu] u(p),
  /// P = p' + p, exact for external nucleons of any masses.
  ///
  /// Index order: 0 = outgoing nucleon (barred), 1 = incoming nucleon.
  /// Parameters:
  ///   BOSON      22, 23 or 24 (default 24 if the legs are p,n, else 22)
  ///   VECTOR_FF  0, 1 (pointlike), 10 (Kelly, default)
  ///   AXIAL_FF   0, 1 (g_A), 20 (Meyer deuterium, default), 21 (lattice)
  ///   and the FF_N_N parameters.
  /// m_N = (m+m')/2, m_pi (pi+- for W, pi0 for Z), s_W^2 and V_ud = |CKM_0_0|
  /// are Sherpa's own values.
  ///
  /// The getter refuses a nucleon beam that carries a PDF: the current
  /// describes the whole nucleon, so its beam must be PDF_SET: Intact.
  ///
  /// m_vanishes flags a current that is identically zero: photon exchange
  /// on a neutron with VECTOR_FF = 0 or 1, where F1n = F2n = 0. It only
  /// warns, and is returned by Vanishes(). Exchange_ME stops the run if
  /// every exchange of a Currents: block vanishes, so a Z exchange next to
  /// the photon keeps a pointlike-neutron run going.
  //////////////////////////////////////////////////////////////////////////
  class VA_N_N : public Current_Base {
  private:
    FormFactor_Base *p_F1p, *p_F2p, *p_F1n, *p_F2n, *p_FA;
    kf_code m_boson;
    double  m_mN, m_mpi, m_s2, m_T3, m_norm;
    bool    m_isp, m_vanishes;

    FormFactor_Base * MakeFF(const std::string & name,const int ffmodel,
                             GeneralModel & model);
  public:
    VA_N_N(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_N_N();
    void CheckTopology() const override;
    inline bool Vanishes() const override { return m_vanishes; }
    void SetModelParameters(struct GeneralModel model) override;
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti) override;
  };
}

#endif
