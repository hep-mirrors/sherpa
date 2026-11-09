#ifndef METOOLS_Current_Library_VA_N_NPi_H
#define METOOLS_Current_Library_VA_N_NPi_H

#include "METOOLS/HadronCurrents/Current_Base.H"

namespace METOOLS {
  class FF_N_N;
  class FF_N_NPi;

  class VA_N_NPi : public Current_Base {
  private:
    FF_N_N   *p_NP[4], *p_CNP[4], *p_F1V;
    FF_N_NPi *p_gpiNN, *p_FPiPole;
    ATOOLS::Flavour_Vector m_ffflavs;
    kf_code m_boson;
    int     m_pipolearg;
    double  m_mN, m_norm, m_s2, m_gA, m_fpi;
    double  m_mNP, m_mCNP, m_mex2, m_mpole2;
    double  m_cNP, m_cCNP, m_cPF, m_cX, m_cCT;
    bool    m_onPF, m_wardcheck;

    double CPi(const int Qin,const int Qout) const;
    void   MakeVertexFF(FF_N_N * ff[4],const int iout,const int iin,
                        const int vff,const int aff,GeneralModel & md);
    void   CheckWard(const ATOOLS::Vec4D & q,const ATOOLS::Vec4C amp[2][2]);
  public:
    VA_N_NPi(const ATOOLS::Flavour_Vector& flavs,
             const std::vector<int>& indices,
             const std::string& name);
    ~VA_N_NPi();
    void CheckTopology() const override;
    void SetModelParameters(struct GeneralModel model) override;
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti) override;
  };
}

//////////////////////////////////////////////////////////////////////////
/// Single-pion production off a nucleon, N(p2) + X(q) -> N'(p4) + pi(p5),
///   X = gamma, Z, W+/-

///
///   NP  : ubar4 V_piNN S_N(p2+q) V^mu_XNN u2      (pion after absorption)
///   CNP : ubar4 V^mu_XNN S_N(p4-q) V_piNN u2      (pion before absorption)
///   PF  : ubar4 V_piNN u2 S_pi(p2-p4) V^mu_Xpipi
///   CT  : ubar4 V^mu_XpiNN u2                     (arXiv:2610.40866, (94))
///
/// with
///      V_piNN = i c_pi g_piNN/(2m_N) (pslash_Nout - pslash_Nin) g5
///         S_N = i(pslash_N + m_N)/(p_N^2-m_N^2)
///        S_pi = i/(p_pi^2-m_pi^2)
///    V^mu_XNN = i Gamma^mu = i (F1 g^mu + F2/(2m_N) i sigma^{mu nu} q_nu - FA g^mu g5 - FP/(2m_N) q^mu g5)
///                  the VA_N_N vertex, q = p_Nout - p_Nin, Q^2 = -q^2, form factors from FF_N_N:
///                    photon : F_i = F_i^N,                     F_A = F_P = 0
///                    Z      : F_i = T3 F_i^V - 2 s_W^2 F_i^N,  F_A -> T3 F_A
///                    W      : F_i = F_i^V, times V_ud
///                  F^V = F^p - F^n, F_P = 4 m_N^2/(Q^2+m_pi^2) F_A (PCAC).
///                  No Gordon identity (as in VA_N_N): one nucleon is internal, off shell.

/// Delta diagrams (DeltaP, CDeltaP) and Resonance diagrams (NprimeP, CNprimeP) are not yet included.
///
/// REDUCTION TO XYZ FUNCTIONS.
///   ubar4 p4slash = m4 ubar4,  p2slash u2 = m2 u2,
///   i sigma^{mu nu} q_nu = q^mu - g^mu qslash = -q^mu + qslash g^mu,
///   p5slash g^mu p5slash = 2 p5^mu p5slash - m_pi^2 g^mu,
///   q = p4 + p5 - p2 (exact for m_p != m_n)

///   Defn: chi = {1, g5} for vector/axial,
///   Y : ubar4 chi u2           X : ubar4 p5slash chi u2
///   L : ubar4 g^mu chi u2      W : ubar4 g^mu p5slash chi u2
///   Y, X, L, W = chi 1;  Y5, X5, L5, W5 = chi g5
///   J^mu = m_norm (J_NP + J_CNP + J_PF + J_CT), per (h0,h1)
///   in Calc: Y  = F.Y(0,h0,1,h1,1.,1.),     Y5 = F.Y(0,h0,1,h1,1.,-1.)
///            X  = F.X(0,h0,p5,1,h1,1.,1.),  X5 = F.X(0,h0,p5,1,h1,1.,-1.)
///            L  = F.L(0,h0,1,h1,1.,1.),     L5 = F.L(0,h0,1,h1,1.,-1.)
///            W  = F.W(0,h0,p5,1,h1,1.,1.),  W5 = F.W(0,h0,p5,1,h1,1.,-1.)
///            amp[h0][h1] = m_norm*J
///
/// THE FOUR CURRENTS (i^3 = -i from the three Feynman rules of each).
///
///  NP : N(p2) + X -> N(p) -> N(p4) + pi(p5),  p = p4+p5, intermediate mass m
///    J_NP = -i G D ubar4 (-p5slash) g5 (pslash + m) Gamma^mu u2
///         = -i G D ubar4 g5 (A + B p5slash) Gamma^mu u2
///    A = p^2 - m4^2,  B = m4 + m,  D = 1/(p^2 - m^2),  G = c_pi g_piNN(-m_pi^2)/(2m_N)
///    term by term, with P(y,l,w) = -A *y + 2B p5^mu *l - B *w
///      F1 g^mu               ->  F1 P(Y5,L5,W5)
///      -FA g^mu g5           -> -FA P(Y,L,W)
///      -FP/(2m_N) q^mu g5    -> -FP/(2m_N) (A Y - B X) q^mu
///      F2/(2m_N) (q^mu - g^mu qslash), with qslash u2 = (p4slash + p5slash - m2) u2
///                            ->  F2/(2m_N) [ (A Y5 - B X5) q^mu + m2 P(Y5,L5,W5)
///                                  - A (2 p4^mu Y5 - m4 L5)
///                                  - B (-2 p4^mu X5 + 2 p4.p5 L5 - m4 (2 p5^mu Y5 - W5))
///                                  - (A W5 - 2B p5^mu X5 + B m_pi^2 L5) ]
///    full current:
///      J_NP = -i G D [ F1 P(Y5,L5,W5) - FA P(Y,L,W) + F2/(2m_N) [F2 bracket]
///                      - FP/(2m_N) (A Y - B X) q^mu ]
///    in Calc (A, B, GD = G D, m_mNP = m, m2N = 2m_N):
///      P   = -A*l + (2.*B*y)*p5 - B*w
///      sig = (A*Y5-B*X5)*q + m2*P(Y5,L5,W5)
///            - ( A*((2.*Y5)*p4-m4*L5) + B*((-2.*X5)*p4+2.*p4p5*L5-m4*((2.*Y5)*p5-W5)) )
///            - ( A*W5-(2.*B*X5)*p5+B*mpi2*L5 )
///      J  += Complex(0.,-GD)*( F1*P(Y5,L5,W5) - FA*P(Y,L,W) + (F2/m2N)*sig - (FP/m2N*(A*Y-B*X))*q )
///
///  CNP : N(p2) -> N(p') + pi(p5), N(p') + X -> N(p4),  p' = p2-p5, intermediate mass m'
///    J_CNP = -i G D' ubar4 Gamma^mu (p'slash + m') (-p5slash) g5 u2
///          = -i G D' ubar4 Gamma^mu (A' - B' p5slash) g5 u2
///    A' = p'^2 - m2^2,  B' = m' + m2,  D' = 1/(p'^2 - m'^2)
///    term by term, with P'(l,w) = A' *l - B' *w:
///      F1' g^mu              ->  F1' P'(L5,W5)
///      -FA' g^mu g5          -> -FA' P'(L,-W)
///      -FP'/(2m_N) q^mu g5   -> -FP'/(2m_N) (A' Y + B' X) q^mu
///      F2'/(2m_N) (-q^mu + qslash g^mu),
///                            ->  F2'/(2m_N) [ -(A' Y5 - B' X5) q^mu + m4 P'(L5,W5)
///                                   + A' (2 p5^mu Y5 - W5) - B' (2 p5^mu X5 - m_pi^2 L5)
///                                   - A' (2 p2^mu Y5 + m2 L5)
///                                   + B' (2 p2^mu X5 - 2 p2.p5 L5 - m2 W5) ]
///    full current:
///      J_CNP = -i G D' [ F1' P'(L5,W5) - FA' P'(L,-W) + F2'/(2m_N) [F2' bracket]
///                        - FP'/(2m_N) (A' Y + B' X) q^mu ]
///    in Calc (Ac = A', Bc = B', GDc = G D', m_mCNP = m', F1c..FPc = F1'..FP'):
///      sig = (-(Ac*Y5-Bc*X5))*q + m4*(Ac*L5-Bc*W5)
///            + Ac*((2.*Y5)*p5-W5) - Bc*((2.*X5)*p5-mpi2*L5)
///            - ( Ac*((2.*Y5)*p2+m2*L5) - Bc*((2.*X5)*p2-2.*p2p5*L5-m2*W5) )
///      J  += Complex(0.,-GDc)*( F1c*(Ac*L5-Bc*W5) - FAc*(Ac*L+Bc*W) + (F2c/m2N)*sig
///                               - (FPc/m2N*(Ac*Y+Bc*X))*q )
///
///  PF : X hits the exchanged pion,  t = (p2-p4)^2
///    J_PF = -i G_PF Y5 (p2 - p4 + p5)^mu,      [ubar4 (p4slash - p2slash) g5 u2 = (m4+m2) Y5]
///    G_PF = c_pi g_piNN(-t)/(2m_N) (m4+m2) c_X F1V / (t - m_pi,ex^2)
///    full current:
///      J_PF = -i c_pi g_piNN(-t)/(2m_N) (m4+m2) c_X F1V / (t - m_pi,ex^2) Y5 (p2 - p4 + p5)^mu
///    in Calc (m_cPF = c_pi, m_cX = c_X, m_mex2 = m_pi,ex^2, K = p2-p4+p5):
///      GPF = m_cPF*p_gpiNN->FF_Assigner(-t)/m2N*(m4+m2)/(t-m_mex2)*m_cX*F1V
///      J  += (Complex(0.,-GPF)*Y5)*K
///
///  CT : contact (Kroll-Ruderman), arXiv:2610.40866 Eq. (94), k_pi = +p5
///    J_CT = -i K_R L5 + i F_R [ L + ((m4-m2) Y + 2 X) q^mu / (2(Q^2 + m_pole^2)) ]
///                                [ubar4 (qslash + p5slash) u2 = (m4-m2) Y + 2 X]
///    K_R = c_hat g_A F1V / f_pi   (x (1 - 2 s_W^2) for Z)
///    F_R = c_hat F_PiPole / f_pi     (Z, W only; F_PiPole at Q^2, or -t for PIPOLE_FF_ARG 1)
///    full current:
///      J_CT = -i c_hat/f_pi [ g_A F1V L5 - F_PiPole ( L + ((m4-m2) Y + 2 X) q^mu / (2(Q^2 + m_pole^2)) ) ]
///      (g_A term x (1 - 2 s_W^2) for Z; F_PiPole term for Z, W only)
///    in Calc (m_cCT = c_hat, m_mpole2 = m_pole^2, m_s2 = s_W^2, p_FPiPole = F_PiPole):
///      KR  = m_cCT*m_gA/m_fpi*F1V*(m_boson==kf_Z ? 1.-2.*m_s2 : 1.)
///      FR  = m_cCT/m_fpi*p_FPiPole->FF_Assigner(m_pipolearg ? -t : Q2)
///      J  += Complex(0.,-KR)*L5
///      J  += Complex(0.,FR)*( L+(((m4-m2)*Y+2.*X)/(2.*(Q2+m_mpole2)))*q )   (not photon)
///
/// ISOSPIN AND CHARGE (algorithmic, from Flavour::Charge()):
///   Q_X = Q4 + Q_pi - Q2; a diagram is kept only if every vertex is
///   charge allowed: NP needs Q2+Q_X in {0,1}, CNP needs Q2-Q_pi in {0,1}.
///   c_pi  = +1 (p, pi0), -1 (n, pi0), sqrt2 (charged pion)
///   c_X   (PF) : photon Q_pi, Z Q_pi (1-2 s_W^2),
///                W sqrt2 c_W, c_W = +1 (W pi^-+ -> pi0), -1 (W pi0 -> pi^+-)
///   c_hat (CT) : photon, Z Q_pi/sqrt2; W 1 for a pi0, else
///                +-1/sqrt2 (W+: n -> n pi+ +, p -> p pi+ -; W- by p<->n).
/// The couplings e, g/(2c_W), g/(2sqrt2) are on the Exchange_ME line;
/// the W carries |V_ud| = |CKM_0_0| here, as in VA_N_N.
///
/// Index order: 0 = outgoing nucleon, 1 = incoming nucleon, 2 = pion.
/// Parameters:
///   BOSON       22, 23, 24 (default 24 if the charges need a W, else 22)
///   VECTOR_FF, AXIAL_FF   as for VA_N_N (FF_N_N models)
///   PION_NN_FF  FF_N_NPi FpiNN model (default 1, Goldberger-Treiman)
///   PIPOLE_FF      FF_N_NPi FPiPole model: 0, 1, 11 rho dominance (default), 12 F1V
///   PIPOLE_FF_ARG  0: x = Q^2 (default), 1: x = -t (HNV)
///   GA (1.2754), fpi (0.093)
///
/// m_wardcheck: if true, the first photon event checks q.J = 0 and
/// prints |q.J|/(|q||J|). Exact for m_p = m_n; with Sherpa's masses the
/// charged channels carry O((m_n-m_p)/m_N) differences.
//////////////////////////////////////////////////////////////////////////

#endif
