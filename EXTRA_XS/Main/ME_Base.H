#ifndef EXTRA_XS_Main_ME_Base_H
#define EXTRA_XS_Main_ME_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Main/Spin_Structure.H"

namespace EXTRAXS {

  class ME_Base {
  protected:
    const PHASIC::Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    METOOLS::Spin_Amplitudes m_res;
    METOOLS::Amplitude_Tensor m_colres;
  public:
    ME_Base(const PHASIC::Process_Info& pi,
            const ATOOLS::Flavour_Vector& flavs, const Particle_Vector& parts);
    virtual ~ME_Base();

    virtual void Calc(const ATOOLS::Vec4D_Vector& momenta) = 0;
    inline const METOOLS::Spin_Amplitudes& Result() const { return m_res; }
    inline const METOOLS::Amplitude_Tensor& Resultcol() const { return m_colres; }

    static ME_Base* GetME(const PHASIC::Process_Info& pi);
    static ME_Base* GetME(const std::string& tag,
                          const PHASIC::Process_Info& pi);
  };
}

#define DECLARE_ME_GETTER(NAME,TAG)                                \
  DECLARE_GETTER(NAME,TAG,ME_Base,Process_Info);                   \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
