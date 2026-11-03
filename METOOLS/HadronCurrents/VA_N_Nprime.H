#ifndef METOOLS_Current_Library_VA_N_Nprime_H
#define METOOLS_Current_Library_VA_N_Nprime_H

#include "METOOLS/HadronCurrents/Current_Base.H"

namespace METOOLS {
//////////////////////////////////////////////////////////////////////////
/// Electromagnetic transition current between on-shell spin-1/2 baryons.
///
///   J^mu = i ubar(p') [ F1prime(Q^2) T1^mu
///                       + F2prime(Q^2) T2^mu ] GammaP u(p)
///
/// where
///
///   T1^mu = gamma^mu - q^mu qslash/q^2
///   T2^mu = i sigma^{mu nu} q_nu/(m_R + m_N)
///
/// q = p' - p and Q^2 = -q^2. The form factors F1prime and
/// F2prime describe the baryon transition and are not, in general,
/// the elastic nucleon form factors.
///
/// GammaP specifies the parity structure of the transition and must
/// be consistent with the parities of the initial and final states.
///
/// Both T1^mu and T2^mu are transverse, so q_mu J^mu = 0,
/// provided the vertex contains no additional non-transverse terms.
///
/// The overall phase and normalisation must be consistent with the
/// coupling convention used by Exchange_ME.
///
/// Index order: 0 = outgoing baryon (barred), 1 = incoming baryon.
///
/// Parameters:
///   BOSON      22 (photon)
///   VECTOR_FF  Model selector for the transition form factors.
///
/// The transition form factors and any required masses or normalisation
/// parameters must be supplied by the resonance form-factor model.
//////////////////////////////////////////////////////////////////////////
  class VA_N_Nprime : public Current_Base {
  private:
    FormFactor_Base *p_F1prime, *p_F2prime;
    kf_code m_boson;
    double  m_N, m_R, m_sum, m_diff, m_norm;
    bool    m_vanishes;

    FormFactor_Base * MakeFFprime(const std::string & name,const int ffmodel,
                             GeneralModel & model);
  public:
    VA_N_Nprime(const ATOOLS::Flavour_Vector& flavs,
           const std::vector<int>& indices,
           const std::string& name);
    ~VA_N_Nprime();
    void CheckTopology() const override;
    inline bool Vanishes() const override { return m_vanishes; }
    void SetModelParameters(struct GeneralModel model) override;
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti) override;
  };
}

#endif
