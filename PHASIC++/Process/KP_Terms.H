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
    ATOOLS::subscheme::code    m_subtype;

    Process_Base    *p_proc;
    PDF::NLOMC_Base *p_nlomc;
    Massive_Kernels *p_kernel;

    MODEL::Coupling_Data *p_cpl;

    ATOOLS::Flavour_Vector m_flavs;

    bool   m_massive, m_cemode, m_negativepdf;
    double m_kpca[8], m_kpcb[8], m_kpcta[8], m_kpctb[8];
    double m_cpldef, m_NC;
    double m_epsa, m_epsb, m_epsdefa, m_epsdefb;
    int    m_Vsubmode, m_facscheme;
    bool   m_sa,m_sb;
    int    m_typea, m_typeb;

    std::vector<size_t> m_plist;

    std::vector<xp_coeff> m_xpa, m_xpb, m_xpta, m_xptb;

    void RegisterDefaults() const;
    void SetColourFactors();
    void SetMassive();

    double g(const double &x, const double &muf2);
    double k(const double &x, const double &muf2);

  public:

    KP_Terms(Process_Base *const proc,const ATOOLS::sbt::subtype st,
             const std::vector<size_t>& partonlist);

    ~KP_Terms();

    void SetCoupling(const MODEL::Coupling_Map *cpls);

    void Calculate(const ATOOLS::Vec4D_Vector &mom,
                   const std::vector<std::vector<double > > &dsij,
                   const double &x0,const double &x1,
                   const double &eta0,const double &eta1,
                   const double &etar0,const double &etar1,
                   const double &wgt);

    void CalculateEpsilon(const ATOOLS::Vec4D_Vector &mom,
                          const std::vector<std::vector<double > > &dsij,
                          const double &eta0,const double &eta1,
                          const double &eps0,const double &eps1);;

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

    inline void SetEpsilonA(const double& eps) { m_epsdefa=eps; }
    inline void SetEpsilonB(const double& eps) { m_epsdefb=eps; }
    inline double EpsilonA() const { return m_epsdefa; }
    inline double EpsilonB() const { return m_epsdefb; }

  };// end of class KP_Terms

  // structure of kp coefficients for beam a and b, kpca and kpcb
  //
  // KP = 1/f [
  //            kpc[0] fq + kpc[1] fqx + kpc[2] fg + kpc[3] fgx 
  //            + ( kpc[4] fq + kpc[5] fqx + kpc[6] fg + kpc[7] fgx )
  //              * log(muF2/muF2ref)
  //
  //            + kpct[0] fq + kpct[1] intfq + kpct[2] fg + kpct[3] intfg
  //            + ( kpct[4] fq + kpct[5] intfq + kpct[6] fg + kpct[7] intfg )
  //              * log(muF2/muF2ref)
  //          ]
  //
  // the first set of terms are needed for every PDF,
  // the second set of terms is needed only for PDFs with an integrable 
  // singularity at x = 1 ( x .. PDF momentum fraction )

};// end of namespace PHASIC

#endif
