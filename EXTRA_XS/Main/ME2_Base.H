#ifndef EXTRA_XS_Main_ME2_Base_H
#define EXTRA_XS_Main_ME2_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace EXTRAXS {

  class ME2_Base {
  protected:
    const PHASIC::Process_Info& m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;

    int **p_colours;
    size_t m_oew, m_oqcd, m_sintt;
  public:
    ME2_Base(const PHASIC::Process_Info& pi,
             const ATOOLS::Flavour_Vector& flavs);
    virtual ~ME2_Base();

    virtual double operator()(const ATOOLS::Vec4D_Vector& momenta) = 0;
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);

    inline int **Colours() const { return p_colours; }
    inline size_t OrderEW() const { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }
    inline size_t SIntType() const { return m_sintt; }

    static ME2_Base* GetME2(const PHASIC::Process_Info& pi);
    static ME2_Base* GetME2(const std::string& tag,
                            const PHASIC::Process_Info& pi);
  };
}

#define DECLARE_ME2_GETTER(NAME,TAG)                               \
  DECLARE_GETTER(NAME,TAG,ME2_Base,Process_Info);                  \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
