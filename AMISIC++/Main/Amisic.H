#ifndef AMISIC_Main_Amisic_H
#define AMISIC_Main_Amisic_H

#include "AMISIC++/Perturbative/Single_Collision_Handler.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Tools/Impact_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/File_IO_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

/*!
  \file Amisic.H
  \brief Declares the class Amisic
*/

namespace ATOOLS {
  class Cluster_Amplitude;
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
    double m_sigmaND_norm;
    double m_b, m_bfac, m_residualE1, m_residualE2, m_pt2, m_Enorm;

    MI_Processes           * p_processes;
    Over_Estimator           m_overestimator;
    Impact_Parameter         m_impact;
    Single_Collision_Handler m_singlecollision;

    bool   m_isMinBias, m_isFirst, m_ana;
    size_t m_nscatters;
    std::map<std::string,ATOOLS::Histogram *> m_histos; 
    
    bool InitParameters();
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
		    PDF::ISR_Handler *const isr);
    
    bool InitMPIs(const double & scale);
    bool VetoEvent(const double & scale=-1.);
    int  InitMinBiasEvent(ATOOLS::Blob_List * blobs);
    void Reset();
    void CleanUp();
    void Test();
    bool VetoScatter(ATOOLS::Blob * blob);
    ATOOLS::Blob * GenerateScatter();
    ATOOLS::Cluster_Amplitude * ClusterConfiguration(ATOOLS::Blob * blob);
    
    
    inline void SetMaxEnergies(const double & E1,const double & E2) {
      m_residualE1 = E1;
      m_residualE2 = E2;
      m_singlecollision.SetResidualX(E1/m_Enorm,E2/m_Enorm);
    }
    inline void SetMaxScale(const double & scale) {
      m_pt2 = ATOOLS::sqr(scale);
      m_singlecollision.SetLastPT2(m_pt2);
    }
    inline void SetB(const double & b=-1.) {
      // Select b and set the enhancement factor
      m_b    = (b<0.)?m_impact.SelectB(m_pt2):b;
      m_bfac = ATOOLS::Max(0.,m_impact(m_b));
    }
    inline void SetMinBias(const bool & isMinBias=true) {
      m_isMinBias = isMinBias;
    }
    inline void SetFirst(const bool & first=true) {
      m_isFirst = first;
    }
    inline void SetMassMode(const int & mode) {
      p_processes->SetMassMode(mode);
    }
    inline int ShiftMasses(ATOOLS::Cluster_Amplitude * ampl) {
        return p_processes->ShiftMasses(ampl);
    }
    inline ATOOLS::Vec4D SelectPositionForScatter(const double & b) const {
      return m_impact.SelectPositionForScatter(b);
    }
    inline ATOOLS::Vec4D  SelectRelativePositionForParton() const {
      return m_impact.SelectRelativePositionForParton();
    }
    inline const double ScaleMin()    const { return (*mipars)("pt_min"); }
    inline const double & ScaleMax()  const { return m_pt2; }
    inline const double & PT2()       const { return m_pt2; }
    inline const bool   & IsMinBias() const { return m_isMinBias; }
    inline const bool   & IsFirst()   const { return m_isFirst; }
    inline const double & B()         const { return m_b; }
    inline double Mass(const ATOOLS::Flavour &fl) const {
      return p_processes->Mass(fl);
    }
  };

}

#endif
