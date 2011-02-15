#ifndef EXTRA_XS_Main_ME2_Base_H
#define EXTRA_XS_Main_ME2_Base_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Coupling_Data.H"

namespace EXTRAXS {

  class ME2_Base: public PHASIC::Tree_ME2_Base {
  protected:
    int **p_colours;
    size_t m_oew, m_oqcd, m_sintt;
    MODEL::Coupling_Data *p_aqcd, *p_aqed;
    double m_symfac;
  public:
    ME2_Base(const PHASIC::Process_Info& pi,
             const ATOOLS::Flavour_Vector& flavs);
    virtual ~ME2_Base();

    void SetCouplings(MODEL::Coupling_Map *const cpls);

    double CouplingFactor(const int oqcd,const int oew) const;

    double Calc(const ATOOLS::Vec4D_Vector &p);

    virtual double operator()(const ATOOLS::Vec4D_Vector& momenta) = 0;
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);

    inline int **Colours() const { return p_colours; }
    inline size_t OrderEW() const { return m_oew;  }
    inline size_t OrderQCD() const { return m_oqcd; }
    inline size_t SIntType() const { return m_sintt; }

  };
}

#endif
