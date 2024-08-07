#ifndef AMISIC_Main_Amisic_H
#define AMISIC_Main_Amisic_H

#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
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

namespace ATOOLS {
  class Cluster_Amplitude;
  class MI_Processes;
}

namespace AMISIC {
  /*!
    \namespace AMISIC
    \brief The module for the generation of a perturbative underlying
    event according to the Sjostrand-van der Zijl model.

    AMISIC++ is an accronym for <em>A</em> <em>M</em>ultiple
    <em>I</em>nteraction <em>S</em>imulation <em>I</em>n
    <em>C++</em>.

    The AMISIC namespace contains the module for the generation of a
    perturbative underlying event according to the Sjostrand-van der Zijl
    model.  The original model is based on perturbative QCD 2->2 scatters
    modified through a simple IR regularisation and convoluted with regular
    PDFs and and a matter-overlap function between the incident hadronic
    states which gives rise to an interaction probability.

    In SHERPA we added/plan to add other 2->2 scatters, for example for photon
    and quarkonia production.
  */

  class Amisic: public ATOOLS::File_IO_Base, ATOOLS::Mass_Selector {
  private:
    MI_Processes             * p_processes;
    Hadronic_XSec_Calculator * p_xsecs;
    Over_Estimator             m_overestimator;
    Impact_Parameter           m_impact;
    Single_Collision_Handler   m_singlecollision;

    mitype m_type;
    bool         m_variable_s;
    double       m_sigmaND_norm;
    double       m_b, m_bfac, m_S, m_Y, m_pt2;

    bool         m_isMinBias, m_isFirst, m_ana;
    size_t       m_nscatters;
    std::map<std::string,ATOOLS::Histogram *> m_histos;

    void InitParametersAndType(PDF::ISR_Handler *const isr);
    void UpdateForNewS();
    void CreateAmplitudeLegs(ATOOLS::Cluster_Amplitude * ampl,
			     ATOOLS::Blob * blob);
    void FillAmplitudeSettings(ATOOLS::Cluster_Amplitude * ampl);
    void InitAnalysis();
    void FinishAnalysis();
    void Analyse(const bool & last);
 public:
    Amisic();
    ~Amisic();

    bool Initialize(MODEL::Model_Base *const model,
		    PDF::ISR_Handler *const isr,
		    YFS::YFS_Handler *const yfs,
                    REMNANTS::Remnant_Handler * remnant_handler);

    bool InitMPIs(PDF::ISR_Handler *const isr, const double & scale);
    bool VetoEvent(const double & scale=-1.) const;
    int  InitMinBiasEvent();
    int  InitRescatterEvent();
    void SetB(const double & b=-1.);

    void Reset();
    void CleanUpMinBias();
    void Test();
    ATOOLS::Blob * GenerateScatter();
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob * blob);

    void Update(const PDF::ISR_Handler *isr);
    inline void SetMaxEnergies(const double & E0,const double & E1) {
      m_singlecollision.SetResidualX(E0,E1);
    }
    inline void SetMaxScale(const double & scale) {
      m_pt2 = ATOOLS::sqr(scale);
      m_singlecollision.SetLastPT2(m_pt2);
    }
    inline void SetMassMode(const int & mode) {
      p_processes->SetMassMode(mode);
    }
    inline int ShiftMasses(ATOOLS::Cluster_Amplitude * ampl) {
        return p_processes->ShiftMasses(ampl);
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
    inline ATOOLS::Vec4D  SelectRelativePositionForParton() const {
      return m_impact.SelectRelativePositionForParton();
    }
    inline const double ScaleMin()    const { return (*mipars)("pt_min"); }
    inline const double ScaleMax()    const { return sqrt(m_pt2); }
    inline const bool   & IsMinBias() const { return m_isMinBias; }
    inline const double & B()         const { return m_b; }
    inline double Mass(const ATOOLS::Flavour &fl) const {
      return p_processes->Mass(fl);
    }
    inline void CleanUp()                   { m_isFirst = true; }
  };

}

#endif
