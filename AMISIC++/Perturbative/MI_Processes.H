#ifndef AMISIC_Perturbative_MI_Processes_H
#define AMISIC_Perturbative_MI_Processes_H

#include "AMISIC++/Perturbative/MI_Integrator.H"
#include "AMISIC++/Perturbative/MI_Process_Group.H"
#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "AMISIC++/Tools/Matter_Overlap.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include <array>

namespace REMNANTS { class Remnant_Base; }
namespace BEAM     { class Beam_Spectra_Handler; }

namespace AMISIC {
  class MI_Processes : public PHASIC::ME_Generator_Base {
  private:
    MI_Integrator                           m_integrator;
    std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> p_remnants;
    std::array<PDF::PDF_Base *, 2>          p_pdf;
    MODEL::One_Running_AlphaS             * p_alphaS;
    MODEL::Running_AlphaQED               * p_alpha;
    Hadronic_XSec_Calculator              * p_xsecs;
    Matter_Overlap                        * p_overlap;

    scale_scheme    m_muR_scheme, m_muF_scheme;
    std::string           m_scale_scheme, m_kfactor_scheme;
    double                m_muR_fac, m_muF_fac, m_muR2, m_muF2;
    double                m_pt0, m_pt02, m_ptmin, m_ptmin2, m_ptmax2, m_Emin;
    std::array<double, 2> m_xmin, m_xmax, m_resx;
    double                m_ecms, m_S, m_S_lab, m_b, m_lastxs;
    double                m_xshard;
    bool                  m_dynamic;
    std::list<MI_Process_Group * >   m_groups;

    bool   InitializeAllProcesses();
    void   SetPDFs();
    void   SetAlphaS();
    void   CalcPDFs(const double & x1,const double & x2);

    MI_Process * SelectProcess();
    void   Test();
  public:
    explicit MI_Processes();
    ~MI_Processes();

    bool   Initialize(MODEL::Model_Base *const model,
		      BEAM::Beam_Spectra_Handler *const beam,
		      PDF::ISR_Handler *const isr,
            YFS::YFS_Handler *const yfs);
    bool   CheckCrossSection();
    void   UpdateS(const double & s);
    /////////////////////////////////////////////////////////////////////////
    // The sum over the PDF-weighted matrix elements, grouped by parton
    // content:
    // [x_1 f_i(x_1, mu^2) x_2 f(x_2, mu^2)]  [pi / shat^2] *
    // {alpha_S, alpha alpha_S, alpha}^2 |M_ij(shat,that,uhat)|^2
    /////////////////////////////////////////////////////////////////////////
    double operator()(const bool & withOverlap=true);
    double operator()(const double & shat,const double & that,
		      const double & uhat,
		      const double & x1,const double & x2,
		      const bool & withOverlap=true);
    void   CalcScales(const double & shat,const double & that,
		      const double & uhat);
    ATOOLS::Blob * FillHardScatterBlob();
    /////////////////////////////////////////////////////////////////////////
    // Return dsigmahat/dpT^2 (in GeV^-4) for a fixed value of pT^2 and, if
    // the matter overlap is dynamic, a fixed value of b.
    // Essentially an integral [dx_1 dx_2 dsigmahat/dpT^2].
    /////////////////////////////////////////////////////////////////////////
    inline double dSigma(const double & pt2) {
      return m_integrator.XSecInPT2Bin(pt2, m_S);
    }
    /////////////////////////////////////////////////////////////////////////
    // Returns the total parton-level cross section sigmahat (in GeV^-2)
    // for a fixed value of s and, if the matter overlap is dynamic, a
    // fixed value of b.
    // Essentially an integral [dpT^2 dx_1 dx_2 dsigmahat/dpT^2].
    /////////////////////////////////////////////////////////////////////////
    inline double SigmaTot(const bool & withMO=true, const double & s=-1.) {
      return m_integrator((s>=0?s:m_S),withMO,false);
    }
    /////////////////////////////////////////////////////////////////////////
    // Ratio of parton-level and non-diffractive cross section.
    /////////////////////////////////////////////////////////////////////////
    inline double XSratio(const double & s) { return p_xsecs->XSratio(s); }

    inline const double  & Ecms()   const { return m_ecms; }
    inline const double  & S()      const { return m_S; }
    inline const double  & PT02()   const { return m_pt02; }
    inline const double  & PTMin()  const { return m_ptmin; }
    inline const double  & PT2Min() const { return m_ptmin2; }
    inline const double  & EMin()   const { return m_Emin; }
    inline const double  & MuR2()   const { return m_muR2; }
    inline const double  & MuF2()   const { return m_muF2; }
    inline const size_t    Size()   const { return m_groups.size(); }
    inline const double  & B()      const { return m_b; }
    inline const double  & Bmax()   const { return p_overlap->Bmax(); }
    inline MI_Integrator * GetIntegrator() { return &m_integrator; }
    inline MODEL::One_Running_AlphaS * AlphaS()     const { return p_alphaS; }
    inline Hadronic_XSec_Calculator  * GetXSecs()   const { return p_xsecs; }inline Matter_Overlap * GetOverlap() const { return p_overlap; }
    inline std::shared_ptr<REMNANTS::Remnant_Base> GetRemnant(const size_t & beam) const {
      return p_remnants[beam];
    }
    inline void SetB(const double & b) { m_b = b; }
    inline void SetMassMode(const int & mode) {
      if (mode==0) abort();
      m_massmode = 1;
    }
    inline void SetXSecCalculator(Hadronic_XSec_Calculator *const xsecs) {
      p_xsecs = xsecs;
    }
    inline void SetMatterOverlap(Matter_Overlap *const mo) {
      p_overlap = mo;
      m_dynamic = p_overlap->IsDynamic();
    }
    inline void SetResidualX(const double & r1,const double & r2) {
      m_resx[0] = r1; m_resx[1] = r2;
    }
    inline int  PerformTests() { return 0; }
    inline bool NewLibraries() { return false; }
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &,
					    bool add) final {
      THROW(fatal_error,
	    "MI_Processes::InitializeProcess should not be called.");
    }
  };
}

#endif
