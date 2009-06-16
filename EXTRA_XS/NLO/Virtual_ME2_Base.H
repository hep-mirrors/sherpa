#ifndef EXTRA_XS_NLO_Virtual_ME2_Base_H
#define EXTRA_XS_NLO_Virtual_ME2_Base_H

#include "EXTRA_XS/Main/ME2_Base.H"
#include "HELICITIES/Loops/Divergence_Array.H"

namespace EXTRAXS {

  class Virtual_ME2_Base {
  protected:
    const PHASIC::Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    HELICITIES::DivArrD m_res;
    bool m_needsborn;
    double m_born;
    double m_mur2;
  public:
    Virtual_ME2_Base(const PHASIC::Process_Info& pi,
                  const ATOOLS::Flavour_Vector& flavs);
    virtual ~Virtual_ME2_Base();

    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta)=0;
    inline const HELICITIES::DivArrD& Result() { return m_res; }

    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);

    virtual double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    inline double ME_Finite() { return m_res.GetFinite(); }
    inline double ME_E1() { return m_res.GetIR(); }
    inline double ME_E2() { return m_res.GetIR2(); }

    inline bool NeedsBorn() const { return m_needsborn; }
    inline void SetBorn(const double& born) { m_born=born; }
    inline void SetRenScale(const double& mur2) { m_mur2=mur2; }

    static Virtual_ME2_Base* GetME2(const PHASIC::Process_Info& pi);
    static Virtual_ME2_Base* GetME2(const std::string& tag,
                                    const PHASIC::Process_Info& pi);
  };

  class Trivial_Virtual : public Virtual_ME2_Base {
  public:
    Trivial_Virtual(const PHASIC::Process_Info& pi,
                 const ATOOLS::Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs) {}

    inline void Calc(const ATOOLS::Vec4D_Vector& momenta) {}
  };
}

#define DECLARE_VIRTUALME2_GETTER(NAME,TAG) \
  DECLARE_GETTER(NAME,TAG,Virtual_ME2_Base,Process_Info); \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
