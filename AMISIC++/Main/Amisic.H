#ifndef AMISIC_Main_Amisic_H
#define AMISIC_Main_Amisic_H

#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
#include "AMISIC++/Perturbative/Sudakov_Argument.H"
#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Tools/Impact_Parameter.H"
#include "AMISIC++/Tools/Hadronic_XSec_Calculator.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "YFS/Main/YFS_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

/*!
  \file Amisic.H
  \brief Declares the class Amisic
*/

namespace ATOOLS { class Cluster_Amplitude; }

namespace AMISIC {
  /*!
    \namespace AMISIC
    \brief The module for the generation of a perturbative underlying event
    according to the Sjostrand-van der Zijl model (Phys.Rev.D 36 (1987) 2019)
    or its extension in the Corke-Sjostrand model (JHEP 05 (2011) 009).

    AMISIC++ is an accronym for <em>A</em> <em>M</em>ultiple
    <em>I</em>nteraction <em>S</em>imulation <em>I</em>n
    <em>C++</em>.

    The AMISIC namespace contains the module for the generation of a
    perturbative underlying event according to the Sjostrand-van der Zijl model
    (Phys.Rev.D 36 (1987) 2019) or its extension in the Corke-Sjostrand model
    (JHEP 05 (2011) 009).
    The original SZ model is based on perturbative QCD 2->2 scatters
    modified through a simple IR regularisation and convoluted with regular
    PDFs and a factorised matter-overlap function between the incident hadronic
    states which gives rise to an interaction probability.  In the CS model
    the matter-overlap depends on the x of the incoming partons initiating the
    scatter, and the available form factors are reduced to a single Gaussian.

    In SHERPA we added/plan to add other 2->2 scatters, for example for photon
    and quarkonia production.
  */

  class Amisic: public ATOOLS::File_IO_Base, ATOOLS::Mass_Selector {
  private:
    Matter_Overlap           m_mo;
    MI_Processes             m_processes;
    Sudakov_Argument         m_sudakovarg;
    Hadronic_XSec_Calculator m_xsecs;
    Over_Estimator           m_overestimator;
    Impact_Parameter         m_impact;
    NonPerturbative_XSecs    * p_soft;
    Single_Collision_Handler m_singlecollision;

    mitype m_type;
    bool         m_variable_s, m_variable_b;
    axis       * p_sbins;
    double       m_sigmaND_norm, m_b, m_bfac, m_S, m_Y, m_pt2;

    bool         m_isMinBias, m_isFirst, m_producedSoft, m_ana;
    double       m_f;
    size_t       m_nscatters;
    unsigned long int m_Nscat, m_Nev;
    std::map<std::string,ATOOLS::Histogram *> m_histos;

    void InitParametersAndType(PDF::ISR_Handler *const isr,
			       REMNANTS::Remnant_Handler * remnant_handler);
    void UpdateForNewS();
    void SetFirstB(const double & b=-1.);
    void SetB(const double & x1=-1.,const double & x2=-1.,
	      const double & scale=-1.);
    void CreateAmplitudeLegs(ATOOLS::Cluster_Amplitude * ampl,
			     ATOOLS::Blob * blob);
    void FillAmplitudeSettings(ATOOLS::Cluster_Amplitude * ampl);
    void InitAnalysis();
    void FinishAnalysis();
    void AnalysePerturbative(const bool & last);
    void Test();
 public:
    Amisic();
    ~Amisic();

    bool Initialize(MODEL::Model_Base *const model,
		    PDF::ISR_Handler *const isr,
		    YFS::YFS_Handler *const yfs,
                    REMNANTS::Remnant_Handler * remnant_handler);

    bool InitMPIs(const double & ptmax,const double & x1,const double & x2,
		  const double & scale);
    bool VetoEvent(const double & scale=-1.) const;
    int  InitMinBiasEvent();
    int  InitRescatterEvent();

    void Reset();
    void CleanUpMinBias();
    void Update(const PDF::ISR_Handler *isr);
    ATOOLS::Blob * GenerateScatter(const bool & force=false);
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob * blob);

    inline void SetMaxEnergies(const double & E0,const double & E1) {
      m_singlecollision.SetResidualX(E0,E1);
    }
    inline void SetMaxScale(const double & scale) {
      m_pt2 = ATOOLS::sqr(scale);
      m_singlecollision.SetLastPT2(m_pt2);
    }
    inline void SetMassMode(const int & mode) {
      m_processes.SetMassMode(mode);
    }
    inline int ShiftMasses(ATOOLS::Cluster_Amplitude * ampl) {
      if (p_processes) return m_processes.ShiftMasses(ampl);
      return 1;
    }
    inline void SetPosition(const size_t & beam,const ATOOLS::Vec4D & pos) {
      m_singlecollision.SetPosition(beam,pos);
    }
    inline void CalculateImpact() {
      m_b    = m_singlecollision.B();
      m_bfac = ATOOLS::Max(0.,m_impact(m_S, m_b));
    }
    inline ATOOLS::Vec4D SelectPositionForScatter(const double & b) const {
      return m_impact.SelectPositionForScatter(b);
    }
    inline const double ScaleMin()    const { return (*mipars)("pt_min"); }
    inline const double ScaleMax()    const { return sqrt(m_pt2); }
    inline const bool   & IsMinBias() const { return m_isMinBias; }
    inline const bool   & IsSoft()    const { return m_producedSoft; }
    inline const double & B()         const { return m_b; }
    inline double Mass(const ATOOLS::Flavour &fl) const {
      if (p_processes) return m_processes.Mass(fl);
      return fl.Mass();
    }
    inline void CleanUp() { m_isFirst = true; }
  };

}

#endif
