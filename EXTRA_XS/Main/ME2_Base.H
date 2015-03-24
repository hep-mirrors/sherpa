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
    double m_symfac;
    std::map<size_t,ATOOLS::Flavour_Vector> m_cfls;
  public:
    ME2_Base(const PHASIC::Process_Info& pi,
             const ATOOLS::Flavour_Vector& flavs);
    virtual ~ME2_Base();

    double CouplingFactor(const int oqcd,const int oew) const;

    double Calc(const ATOOLS::Vec4D_Vector &p);

    virtual double operator()(const ATOOLS::Vec4D_Vector& momenta) = 0;
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);

    int OrderQCD(const int &id=-1);
    int OrderEW(const int &id=-1);

    inline int **Colours() const { return p_colours; }
    inline size_t SIntType() const { return m_sintt; }

    inline const ATOOLS::Flavour_Vector &
    CombinedFlavour(const size_t &idij) { return m_cfls[idij]; }

  };
}

#endif
