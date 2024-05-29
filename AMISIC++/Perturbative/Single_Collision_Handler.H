#ifndef AMISIC_Perturbative_Single_Collision_Handler_H
#define AMISIC_Perturbative_Single_Collision_Handler_H

#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Over_Estimator.H"
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
    MI_Processes                          * p_processes;
    MI_Integrator                         * p_integrator;
    Hadronic_XSec_Calculator              * p_xsecs;
    Over_Estimator                        * p_overestimator;
    NonPerturbative_XSecs                 * p_soft;

    MI_Process      * p_proc;
    double            m_pt2, m_pt2min;
    double            m_S, m_lastpt2, m_residualx1, m_residualx2, m_Ycms;
    double            m_xt, m_ymax, m_yvol, m_y3, m_y4;
    double            m_shat, m_that, m_uhat, m_x1, m_x2, m_muf2, m_mur2;
    std::list<ATOOLS::Blob *> m_prefabs;

    bool m_ana;
    std::map<std::string, ATOOLS::Histogram*> m_histos;
    void InitAnalysis();
    void FinishAnalysis();
    void AnalyseWeight(const double & weight);
    void Analyse(const double & pt2,ATOOLS::Blob * blob);

    bool SelectPT2();
  public:
    Single_Collision_Handler();
    ~Single_Collision_Handler();

    void Init(
	      MI_Processes * processes, Over_Estimator * overestimator);
    void Init(REMNANTS::Remnant_Handler * remnant_handler,
	      NonPerturbative_XSecs * soft);
    void UpdateSandY(double s, double y);
    void PrefabricateBlob();
    ATOOLS::Blob * NextScatter(const bool & force);

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
    inline void AddPrefabricatedBlob(ATOOLS::Blob * blob) {
      m_prefabs.push_back(blob);
    }
    inline ATOOLS::Vec4D Position(const size_t & beam) const {
      return p_remnants[beam]->Position();
    }
    inline ATOOLS::Vec4D InMomentum(const size_t & beam) const {
      return p_remnants[beam]->InMomentum();
    }
    inline const double B() const {
      return sqrt(ATOOLS::dabs((p_remnants[0]->Position()-
				p_remnants[1]->Position()).Abs2()));
    }
    inline const double & muF2()    const { return m_muf2; }
    inline const double&  muR2() const { return m_mur2; }
    inline const double & LastPT2() const { return m_lastpt2; }
    inline const double & PT2()     const { return m_pt2; }
    inline const double & Y3()      const { return m_y3; }
    inline const double & Y4()      const { return m_y4; }

    void Test(const double & Q2,const long int & n);
  };
}

#endif
