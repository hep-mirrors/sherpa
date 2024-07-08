#ifndef AMISIC_Perturbative_MI_Processes_H
#define AMISIC_Perturbative_MI_Processes_H

#include "AMISIC++/Perturbative/MI_Integrator.H"
#include "AMISIC++/Perturbative/MI_Process_Group.H"
#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
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

    scale_scheme    m_muR_scheme, m_muF_scheme;
    std::string           m_scale_scheme, m_kfactor_scheme;
    double                m_muR_fac, m_muF_fac, m_muR2, m_muF2;
    double                m_pt0, m_pt02, m_ptmin, m_ptmin2, m_ptmax2, m_Emin;
    std::array<double, 2> m_xmin, m_xmax, m_resx;
    double                m_ecms, m_S, m_S_lab, m_lastxs, m_fm2GeV2;
    double                m_xshard;
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
    double TotalCrossSection(const double & s,const bool & output=false);
    void   UpdateS(const double & s);
    /////////////////////////////////////////////////////////////////////////
    // The sum over the PDF-weighted matrix elements, grouped by parton
    // content:
    // [x_1 f_i(x_1, mu^2) x_2 f(x_2, mu^2)]  [pi / shat^2] *
    // {alpha_S, alpha alpha_S, alpha}^2 |M_ij(shat,that,uhat)|^2
    /////////////////////////////////////////////////////////////////////////
    double operator()();
    double operator()(const double & shat,const double & that,
		      const double & uhat,
		      const double & x1,const double & x2);
    void   CalcScales(const double & shat,const double & that,
		      const double & uhat);
    int    FillHardScatterBlob(ATOOLS::Blob *& blob,const double & pt2veto=-1.);
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
    inline MI_Integrator * GetIntegrator() { return &m_integrator; }
    inline MODEL::One_Running_AlphaS * AlphaS()     const { return p_alphaS; }
    inline Hadronic_XSec_Calculator  * GetXSecs()   const { return p_xsecs; }
    inline std::shared_ptr<REMNANTS::Remnant_Base> GetRemnant(const size_t & beam) const {
      return p_remnants[beam];
    }
    inline void SetMassMode(const int & mode) {
      if (mode==0) abort();
      m_massmode = 1;
    }
    inline void SetXSecCalculator(Hadronic_XSec_Calculator *const xsecs) {
      p_xsecs = xsecs;
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
