#ifndef PHASIC_Process_Virtual_ME2_Base_H
#define PHASIC_Process_Virtual_ME2_Base_H

#include "PHASIC++/Process/Process_Info.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "MODEL/Main/Coupling_Data.H"

namespace PHASIC {

  class Virtual_ME2_Base {
  protected:
    const Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    METOOLS::DivArrD m_res;
    bool m_needsborn,m_newlibs;
    double m_born;
    double m_mur2;
    MODEL::Coupling_Data *p_aqcd, *p_aqed;
  public:
    Virtual_ME2_Base(const PHASIC::Process_Info& pi,
                  const ATOOLS::Flavour_Vector& flavs);
    virtual ~Virtual_ME2_Base();

    void SetCouplings(MODEL::Coupling_Map *const cpls);

    double CouplingFactor(const int oqcd,const int oew) const;

    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta)=0;
    inline const METOOLS::DivArrD& Result() { return m_res; }

    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);

    virtual double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    inline double ME_Finite() { return m_res.GetFinite(); }
    inline double ME_E1() { return m_res.GetIR(); }
    inline double ME_E2() { return m_res.GetIR2(); }

    inline bool NeedsBorn() const { return m_needsborn; }
    inline void SetBorn(const double& born) { m_born=born; }
    inline void SetRenScale(const double& mur2) { m_mur2=mur2; }
    inline bool NewLibs() const { return m_newlibs; }

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
