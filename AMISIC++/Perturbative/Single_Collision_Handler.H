#ifndef AMISIC_Perturbative_Single_Collision_Handler_H
#define AMISIC_Perturbative_Single_Collision_Handler_H

#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Over_Estimator.H"
#include "AMISIC++/Tools/Interaction_Probability.H"
#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/NonPerturbative_XSecs.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace AMISIC {
  class Single_Collision_Handler {
  private:
    evt_type::code m_evttype;
    std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> p_remnants;
    MI_Processes             * p_processes;
    MI_Integrator            * p_integrator;
    Hadronic_XSec_Calculator * p_xsecs;
    NonPerturbative_XSecs    * p_soft;
    Over_Estimator           * p_overestimator;
    Interaction_Probability  * p_pint;
    Matter_Overlap           * p_overlap;

    MI_Process   * p_proc;
    double         m_pt2, m_pt2min, m_b, m_bfac, m_fb;
    double         m_S, m_lastpt2, m_nextpt2, m_residualx1, m_residualx2, m_Ycms;
    bool           m_done;
    ATOOLS::Blob * m_nextblob;
    std::map<double, ATOOLS::Blob *> m_prefabs;

    double MakeHardScatterBlob(ATOOLS::Blob * blob);
    void   DefinePossibleNext();
    int    SelectPT2();
    void   CopyPrefabricatedBlob(ATOOLS::Blob * blob);

    bool m_ana;
    std::map<std::string, ATOOLS::Histogram*> m_histos;
    void InitAnalysis();
    void FinishAnalysis();
    void AnalyseWeight(const double & weight);
    void Analyse(const double & pt2,ATOOLS::Blob * blob);

  public:
    Single_Collision_Handler();
    ~Single_Collision_Handler();

    void Init(MI_Processes * processes, Over_Estimator * overestimator,
	      Interaction_Probability * pint,Matter_Overlap * overlap);
    void Init(REMNANTS::Remnant_Handler * remnant_handler,
	      NonPerturbative_XSecs * soft);
    void UpdateSandY(double s, double y);
    void Reset();
    bool FirstMinBiasScatter(ATOOLS::Blob * first);
    bool FirstMPI(ATOOLS::Blob * signal);
    bool NextScatter(ATOOLS::Blob * blob);

    inline void SetResidualX(const double & E1, const double & E2) {
      m_residualx1 = E1/InMomentum(0)[0];
      m_residualx2 = E2/InMomentum(1)[0];
      p_processes->SetResidualX(m_residualx1,m_residualx2);
    }
    inline void SetLastPT2(const double pt2=-1.) {
      m_lastpt2 = pt2 < 0. ? m_S / 4. : pt2;
    }
    inline void SetPosition(const size_t & beam,const ATOOLS::Vec4D & pos) {
      p_remnants[beam]->SetPosition(pos);
    }
    inline void SetDone(const bool & done) { m_done = done; }
    inline void AddPrefabricatedBlob(const double & pt2,ATOOLS::Blob * blob) {
      m_prefabs[pt2] = blob;
    }
    inline ATOOLS::Vec4D Position(const size_t & beam) const {
      return p_remnants[beam]->Position();
    }
    inline ATOOLS::Vec4D InMomentum(const size_t & beam) const {
      return p_remnants[beam]->InMomentum();
    }
    inline const double & B()       const { return m_b; }
    inline const double & BFac()    const { return m_bfac; }
    inline const double & LastPT2() const { return m_lastpt2; }
    inline const double & PT2()     const { return m_pt2; }
    inline const bool   & Done()    const { return m_done; }

    void Test(const double & Q2,const long int & n);
  };
}

#endif
