#ifndef AMISIC_Perturbative_MI_Processes_H
#define AMISIC_Perturbative_MI_Processes_H

#include "AMISIC++/Perturbative/MI_Process_Group.H"
#include "AMISIC++/Perturbative/Sudakov_Argument.H"
#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace PDF    { class Remnant_Base; }
namespace BEAM   { class Beam_Spectra_Handler; }

namespace AMISIC {
  class MI_Processes : public PHASIC::ME_Generator_Base {
  private:
    PDF::ISR_Handler             * p_isr;
    std::array<PDF::PDF_Base *, 2> p_pdf;
    MODEL::One_Running_AlphaS    * p_alphaS;
    MODEL::Running_AlphaQED      * p_alpha;
    Hadronic_XSec_Calculator     * p_xsecs;
    Sudakov_Argument             * p_sudakov;

    double      m_muFfac, m_pt0, m_pt02, m_ptmin, m_ptmin2, m_ptmax2;
    std::array<double, 2> m_xmin;
    double      m_ecms, m_S, m_S_lab, m_lastxs;
    std::string m_scale_scheme, m_kfactor_scheme;
    size_t      m_pt2bins, m_sbins, m_MCpoints;
    std::list<MI_Process_Group * >   m_groups;

    bool m_variable_s;

    bool InitializeAllProcesses();
    void SetPDFs();
    void SetAlphaS();
    void CalcPDFs(const double & x1,const double & x2,const double & scale);
    bool PrepareSudakovFactor();
    void Test();
  public:
    explicit MI_Processes(const bool & variable_s);
    ~MI_Processes();

    bool Initialize(MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);
    double operator()(const double & shat,const double & that,
			    const double & uhat,const double & x1,const double & x2);
    double dSigma(const double & pt2);
    MI_Process * SelectProcess();
    void         UpdateS(const double & s);

    inline const double  & Ecms()      const { return m_ecms; }
    inline const double  & S()         const { return m_S; }
    inline const double  & PT02()      const { return m_pt02; }
    inline const double  & PTMin()     const { return m_ptmin; }
    inline const double  & PT2Min()    const { return m_ptmin2; }
    inline const size_t    Size()      const { return m_groups.size(); }
    inline const bool    & VariableS() const { return m_variable_s; }
    inline const double  & PDFXmin(const size_t beam) const { return m_xmin[beam]; }
    inline double XSratio(const double & s)              { return p_sudakov->XSratio(s); }
    inline PDF::PDF_Base * PDF(const size_t beam) const  { return p_pdf[beam]; }
    inline Sudakov_Argument         * GetSudakov() const { return p_sudakov; }
    inline MODEL::One_Running_AlphaS    * AlphaS() const     { return p_alphaS; }
    inline Hadronic_XSec_Calculator * GetXSecs() const   { return p_xsecs; }
    inline void SetMassMode(const int & mode) {
      if (mode==0) abort();
      m_massmode = 1;
    }
    inline void SetXSecCalculator(AMISIC::Hadronic_XSec_Calculator *const xsecs) {
      p_xsecs = xsecs;
    }
    inline int  PerformTests() { return 0; }
    inline bool NewLibraries() { return false; }

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &,
					    bool add) final {
      THROW(fatal_error,"MI_Processes::InitializeProcess should not be called.");
    }
  };
}

#endif
