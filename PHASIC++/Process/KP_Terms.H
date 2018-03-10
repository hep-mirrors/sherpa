#ifndef PHASIC__Process__KP_Terms_H
#define PHASIC__Process__KP_Terms_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Massive_Kernels.H"
#include "MODEL/Main/Coupling_Data.H"
#include "PDF/Main/NLOMC_Base.H"

namespace PDF { class PDF_Base; }

namespace PHASIC {

  struct xp_coeff {
    double xp;
    double kpc;
  };

  class KP_Terms {
  private:
    ATOOLS::sbt::subtype       m_stype;
    ATOOLS::cs_itype::type     m_itype;
    ATOOLS::cs_kcontrib::type  m_kcontrib;

    Process_Base    *p_proc;
    PDF::NLOMC_Base *p_nlomc;
    Massive_Kernels *p_kernel;

    MODEL::Coupling_Data *p_cpl;

    ATOOLS::Flavour_Vector m_flavs;

    bool   m_massive, m_cemode, m_negativepdf;
    double m_kpca[8], m_kpcb[8], m_cpldef, m_NC;
    int    m_Vsubmode, m_collVFF, m_facscheme;
    bool   m_sa,m_sb;
    int    m_typea, m_typeb, m_subtype;

    std::vector<size_t> m_plist;

    std::vector<xp_coeff> m_xpa, m_xpb;

    void SetColourFactors();
    void SetMassive();

  public:

    KP_Terms(Process_Base *const proc,const ATOOLS::sbt::subtype st,
             const std::vector<size_t>& partonlist);

    ~KP_Terms();

    void SetCoupling(const MODEL::Coupling_Map *cpls);

    void Calculate(const ATOOLS::Vec4D_Vector &mom,
                   const std::vector<std::vector<double > > &dsij,
                   const double &x0,const double &x1,
                   const double &eta0,const double &eta1,
                   const double &wgt);

    double Get(PDF::PDF_Base *pdfa, PDF::PDF_Base *pdfb,
               const double &x0,const double &x1,
               const double &eta0,const double &eta1,
               const double &muf02,const double &muf12,
               const double &muf02fac,const double &muf12fac,
               const ATOOLS::Flavour &fl0,const ATOOLS::Flavour &fl1);

    void SetAlpha(const double &aff,const double &afi,
                  const double &aif,const double &aii);
    void SetKappa(const double &kappa);
    void SetNLOMC(PDF::NLOMC_Base *const nlomc);

    void FillMEwgts(ATOOLS::ME_Weight_Info &wgt);

    inline double Coupling() const { return m_cpldef*p_cpl->Factor(); }

    inline Massive_Kernels *Kernel() const { return p_kernel; }

    inline void SetIType(ATOOLS::cs_itype::type itype) { m_itype=itype; }

  };// end of class KP_Terms

};// end of namespace PHASIC

#endif
