#ifndef AMISIC_Perturbative_MI_Processes_H
#define AMISIC_Perturbative_MI_Processes_H

#include "AMISIC++/Perturbative/MI_Process_Group.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace PDF    { class Remnant_Base; }
namespace BEAM   { class Beam_Spectra_Handler; }

namespace AMISIC {
  class MI_Processes : public PHASIC::ME_Generator_Base {
  private:
    PDF::ISR_Handler           * p_isr;
    PDF::PDF_Base              * p_pdf[2];
    MODEL::Model_Base          * p_model;
    MODEL::Running_AlphaS      * p_alphaS;
    MODEL::Running_AlphaQED    * p_alpha;

    double      m_muFfac, m_xmin[2];
    double      m_ecms, m_S, m_pt0, m_pt02, m_ptmin, m_ptmin2, m_ptmax2;
    double      m_xTmin, m_xTmax, m_xTstep;
    double      m_sigmaND;
    double      m_lastxs, m_integral, m_pt2step;
    std::string m_scale_scheme, m_kfactor_scheme;
    int         m_photons, m_quarkonia;
    int         m_nbins, m_MCpoints;
    std::vector<double>             m_diffbins;
    std::vector<double>             m_intbins;
    std::list<MI_Process_Group * >  m_groups;

    bool m_test;

    bool   InitializeAllProcesses();
    void   SetPDFs(); 
    void   SetAlphaS();
    bool   PrepareSudakovFactor();
    double dSigma(const double & pt2);
    void   Test();
    const double SudakovDiffArgument(const double & pt2) const;
  public:
    MI_Processes();
    ~MI_Processes();
    
    bool Initialize(MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    void CalcPDFs(const double & x1,const double & x2,const double & scale);
    const double XSec(const ATOOLS::Vec4D_Vector & momenta);
    const double XSec(const double & shat,const double & that,
		      const double & uhat);
    const double operator()(const double & shat,const double & that,
			    const double & uhat);
    const double SudakovArgument(const double & pt2) const;

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &,
					    bool add) final {
      THROW(fatal_error,
	    "MI_Processes::InitializeProcess should not be called.");
    }
    MI_Process * SelectProcess(); 
  
    inline const double &  Ecms()      const { return m_ecms; }
    inline const double &  PT0()       const { return m_pt0; }
    inline const double &  PTmin()     const { return m_ptmin; }
    inline const double &  XShard()    const { return m_integral; }
    inline const size_t    Size()      const { return m_groups.size(); }
    inline const int    &  MassMode()  const { return m_massmode; }
    inline const double &  PDFXmin(const size_t beam) const { return m_xmin[beam]; }
    inline PDF::PDF_Base * PDF(const size_t beam)  { return p_pdf[beam]; }
    inline MODEL::Running_AlphaS * AlphaS() const  { return p_alphaS; }
    inline void SetSigmaND(const double & sigmaND) { m_sigmaND = sigmaND; }
    inline void SetMassMode(const int & mode) {
      if (mode==0) abort();
      m_massmode = 1;
    }
    inline int  PerformTests() { return 0; }
    inline bool NewLibraries() { return false; }
  };
}

#endif
