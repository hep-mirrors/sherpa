#ifndef METOOLS_Current_Library_VA_P_X_Base_H
#define METOOLS_Current_Library_VA_P_X_Base_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "METOOLS/HadronCurrents/FormFactors/FF_P_X.H"
#include <functional>
#include <map>

namespace METOOLS {

  //////////////////////////////////////////////////////////////////////////
  //
  // Common base for the semileptonic charm transition currents
  //
  //     H^mu = <X(p_X,lambda)| qbar gamma^mu (1-gamma_5) c |D(p)>
  //
  // of charm_semileptonic.tex. One derived class per J^P of X:
  //
  //     VA_P_P   X = 0^-   (K, pi, eta, eta')
  //     VA_P_S   X = 0^+   (K*_0(1430), f0(980), a0(980), f0(500))
  //     VA_P_V   X = 1^-   (K*(892), rho, omega, phi)
  //     VA_P_A   X = 1^+   (K_1(1270), K_1(1400))
  //     VA_P_T   X = 2^+   (K*_2(1430))
  //
  // Internal index convention, matching the existing D-meson entries
  // in Decaydata.yaml:  p_i[0] = decaying pseudoscalar D, p_i[1] = X.
  // The leptonic side is a separate VA_F_F current and the G_F/sqrt(2)
  // prefactor is supplied by Current_ME - so what is built here is
  // exactly V_cq H^mu, nothing more.
  //
  // The resonance X is produced ON its generated invariant mass and is
  // handed to Sherpa's decay handler as a normal unstable particle.
  // Its strong decay and line shape are therefore NOT part of this
  // current (note, "Numerical implementation notes" item 5) - the
  // weak vertex only ever sees the polarization object of X.
  //
  //////////////////////////////////////////////////////////////////////////

  class VA_P_X_Base : public Current_Base {
  protected:
    /// Nominal parent mass M and its square.
    double m_M, m_M2;
    /// Nominal daughter mass (used unless m_dynamicmass is set).
    double m_mX;

    /// V_cq times any isospin/eta-mixing Clebsch and the extra
    /// 1/sqrt(2) for a K_S/K_L final state. Everything that multiplies
    /// H^mu as a pure number lives here.
    double m_norm;
    double m_Vcq, m_isospin;
    bool   m_isKSKL;
    /// If true the stored F(0) already contains |V_cq| (the
    /// "CKM-dressed" f-hat of the note's Sec. on CKM dressing) and
    /// m_Vcq is set to 1 so it is not applied twice.
    bool   m_ckmdressed;

    /// Sign of the Levi-Civita convention relative to
    /// epsilon^{0123}=+1. See the Levi-Civita note in FF_P_X.H: this
    /// is a RELATIVE vector/axial sign and is physical. Default +1,
    /// which maps onto the legacy HADRONS++ charm convention
    /// (fac_costhl = -1 there, with a global i factored differently);
    /// -1 would be the bottom convention.
    double m_epssign;

    /// Use the event-by-event sqrt(s_X) rather than the nominal
    /// resonance mass in the kinematic coefficients (note items 3/4).
    bool   m_dynamicmass;
    /// Floor below which the 1/q^2 structures are replaced by their
    /// analytic limits (note item 1).
    double m_q2min;

    /// Number of recoiling hadrons, p_i.size()-1. 1 for every
    /// quasi-two-body current, 2 for VA_P_PP.
    size_t m_nX;

    /// Assert the expected number of recoiling hadrons. Called by each
    /// derived current so that the message names the actual class.
    void CheckArity(const size_t & nX) const;

    std::map<std::string,FF_P_X*> m_ffs;

    /// Build and register one invariant form factor under 'name'; all
    /// of its parameters are read from the same GeneralModel using the
    /// "<name>_..." key prefix.
    FF_P_X * MakeFF(const std::string & name,GeneralModel & model);
    inline FF_P_X * FF(const std::string & name) { return m_ffs[name]; }

    /// Common YAML keys: Vcq, ISOSPIN, CKM_DRESSED, EPS_SIGN,
    /// DYNAMIC_MASS, Q2MIN.
    void   ReadCommonParameters(GeneralModel & model);

    /// |Delta S| bookkeeping -> default CKM element. Overridable by
    /// the "Vcq" key; the choice is logged at initialisation.
    static int  MesonStrangeness(const ATOOLS::Flavour & fl);
    double      DefaultVcq() const;

    /// N(q^2)/q^2 for a numerator that vanishes at q^2=0 by the
    /// kinematic constraint. Below m_q2min the analytic limit N'(0) is
    /// used, evaluated as a symmetric difference - so the constraint
    /// is respected even for shapes whose derivative is not available
    /// in closed form (note item 1).
    double RegularisedOverQ2(const std::function<double(double)> & N,
                             const double & q2) const;

    /// Three massive spin-1 polarization vectors of an OUTGOING
    /// particle, i.e. already complex conjugated. Taken from
    /// METOOLS::Polarization_Vector so that the basis - and hence the
    /// spin density matrix handed on to the strong decay of X - is the
    /// same one the rest of METOOLS uses.
    void PolarizationVectors(const ATOOLS::Vec4D & p,
                             std::vector<ATOOLS::Vec4C> & eps) const;

    /// Spin-2 polarization tensor contracted with the parent momentum,
    /// eps_T^mu = eps*^{mu nu} p_nu / M  (Eq.effectiveT). Returned for
    /// all five helicities. Built from the spin-1 basis above with the
    /// standard 1 (x) 1 -> 2 Clebsch-Gordan coefficients, so that it
    /// inherits exactly the METOOLS spin-1 phase conventions.
    void EffectiveTensorPolarizations(const ATOOLS::Vec4D & pX,
                                      const ATOOLS::Vec4D & pD,
                                      std::vector<ATOOLS::Vec4C> & epsT) const;

    //////////////////////////////////////////////////////////////////////
    // Shared Lorentz kernels.
    //
    // Each returns H^mu for one Lorentz structure, WITHOUT m_norm and
    // without any global phase, given the parent momentum p and the
    // recoiling hadronic momentum Q (Q^2 = s is always read off the
    // momenta, so a resonance produced on a generated mass and a
    // genuine multi-hadron system are treated identically).
    //
    // Factoring these out is what lets one coherent multi-hadron
    // current reuse the SAME, already validated, decompositions that
    // the quasi-two-body currents use: for a two-pseudoscalar system
    // the P wave is exactly VectorKernel with the vector's
    // polarization replaced by the normalised transverse relative
    // momentum, and the D wave is VectorKernel again with the
    // spin-2 object built from it. Nothing is rederived.
    //////////////////////////////////////////////////////////////////////

    /// 0^- / 0^+ structure: f_+ [P - (Delta/q^2) q] + f_0 (Delta/q^2) q.
    ATOOLS::Vec4C ScalarKernel(const ATOOLS::Vec4D & p,
                               const ATOOLS::Vec4D & Q,
                               FF_P_X * fp,FF_P_X * f0) const;

    /// 1^- / 2^+ structure, Eqs.(PVvector),(PVaxial). 'pol' is the
    /// already-conjugated polarization-like four-vector: eps* for a
    /// vector, eps_T for a tensor, or a wave projector for a
    /// multi-hadron system. A_0(0) is fixed per call by
    /// Eq.(A0constraint) using the generated sqrt(s).
    ATOOLS::Vec4C VectorKernel(const ATOOLS::Vec4D & p,
                               const ATOOLS::Vec4D & Q,
                               const ATOOLS::Vec4C & pol,
                               FF_P_X * V,FF_P_X * A0,
                               FF_P_X * A1,FF_P_X * A2,
                               const double & csign) const;

    /// 1^+ structure, Eqs.(PAvector),(PAaxial): the parity flip of the
    /// above, with (M+m) and (M-m) interchanged.
    ATOOLS::Vec4C AxialKernel(const ATOOLS::Vec4D & p,
                              const ATOOLS::Vec4D & Q,
                              const ATOOLS::Vec4C & pol,
                              FF_P_X * A,FF_P_X * V0,
                              FF_P_X * V1,FF_P_X * V2,
                              const double & csign,
                              const double & phase) const;

    void PrintCommonParameters() const;

  public:
    VA_P_X_Base(const ATOOLS::Flavour_Vector & flavs,
                const std::vector<int> & indices,
                const std::string & name);
    virtual ~VA_P_X_Base();
  };
}

#endif
