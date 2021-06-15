#ifndef EXTRA_XS_Main_ME2_Base_H
#define EXTRA_XS_Main_ME2_Base_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "MODEL/Main/Coupling_Data.H"

namespace EXTRAXS {

  class ME2_Base: public PHASIC::Tree_ME2_Base {
  protected:
    std::vector<std::vector<int> > m_colours;
    size_t m_oew, m_oqcd, m_sintt;
    double m_symfac, m_sprimemin, m_sprimemax;
    bool   m_hasinternalscale;

    double m_internalscale;
    std::map<size_t,ATOOLS::Flavour_Vector> m_cfls;

    std::vector<ATOOLS::Flavour> m_flavours;
    std::vector<ATOOLS::Vec4D>   m_momenta;
  public:
    ME2_Base(const PHASIC::External_ME_Args&);
    virtual ~ME2_Base();

    double CouplingFactor(const int oqcd,const int oew) const;

    double Calc(const ATOOLS::Vec4D_Vector &p);

    virtual double operator()(const ATOOLS::Vec4D_Vector& momenta) = 0;
    virtual bool SetColours(const ATOOLS::Vec4D_Vector& momenta);
    virtual bool FillFinalState(const ATOOLS::Vec4D_Vector &p) { return true; }

    inline size_t NOut() const { return m_flavours.size()-2; }
    inline const ATOOLS::Flavour_Vector & Flavours()    const { return m_flavours; }
    inline const ATOOLS::Vec4D_Vector   & Momenta()     const { return m_momenta; }
    inline const ATOOLS::Flavour_Vector * GetFlavours() const { return &m_flavours; }
    inline const ATOOLS::Vec4D_Vector   * GetMomenta()  const { return &m_momenta; }
    
    int OrderQCD(const int &id=-1);
    int OrderEW(const int &id=-1);

    inline std::vector<std::vector<int> > & Colours() { return m_colours; }
    inline size_t SIntType() const { return m_sintt; }

    inline const ATOOLS::Flavour_Vector &
    CombinedFlavour(const size_t &idij) { return m_cfls[idij]; }

    inline const bool HasInternalScale() const { return m_hasinternalscale; }
    inline const double InternalScale()  const { return m_internalscale; }
    inline const double SPrimeMin() const { return m_sprimemin; }
    inline const double SPrimeMax() const { return m_sprimemax; }
  };
}

#endif
