#ifndef EXTRA_XS_NLO_Loop_ME_Base_H
#define EXTRA_XS_NLO_Loop_ME_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "HELICITIES/Loops/Divergence_Array.H"
#include "HELICITIES/Main/Spin_Structure.H"

namespace EXTRAXS {

  typedef HELICITIES::Spin_Structure<HELICITIES::DivArrC> DivAmps;

  class Loop_ME_Base {
  protected:
    const PHASIC::Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    double m_mur2;
    DivAmps m_res;
  public:
    Loop_ME_Base(const PHASIC::Process_Info& pi,
                 const ATOOLS::Flavour_Vector& flavs);
    virtual ~Loop_ME_Base();

    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta)=0;
    inline const DivAmps& Result() const { return m_res; }

    inline void SetRenScale(const double& mur2) { m_mur2 = mur2; }

    static Loop_ME_Base* GetME(const PHASIC::Process_Info& pi);
    static Loop_ME_Base* GetME(const std::string& tag,
                               const PHASIC::Process_Info& pi);

  };

}

#define DECLARE_LOOPME_GETTER(NAME,TAG) \
  DECLARE_GETTER(NAME,TAG,Loop_ME_Base,Process_Info); \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
