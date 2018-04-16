#ifndef AMISIC_Perturbative_Single_Collision_Handler_H
#define AMISIC_Perturbative_Single_Collision_Handler_H

#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Over_Estimator.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

/*!
  \file Single_Collision_Handler.H
  \brief Declares the class Single_Collision_Handler
*/

namespace AMISIC {
  class Single_Collision_Handler {
  private:
    MI_Processes   * p_processes;
    Over_Estimator * p_overestimator;
    double           m_pt2, m_pt2min;
    ATOOLS::Vec4D    m_pbeam1, m_pbeam2;
    double           m_S, m_lastpt2, m_residualx1, m_residualx2;

    MI_Process     * p_proc;
    double           m_xt, m_ymax, m_yvol, m_y3, m_y4, m_x1, m_x2;
    double           m_shat, m_that, m_uhat;
    double           m_muf2, m_mur2;

    bool m_ana;
    std::map<std::string,ATOOLS::Histogram *> m_histos; 
    void InitAnalysis();
    void FinishAnalysis();
    void AnalyseWeight(const double & weight);
    void Analyse(const double & pt2,ATOOLS::Blob * blob);

    bool SelectPT2(const double & pt2);
    bool SelectRapidities();
    bool CalcXs();
    bool CalcMandelstams();
    ATOOLS::Blob * MakeBlob();
  public:
    Single_Collision_Handler();
    ~Single_Collision_Handler();

    void Init();
    ATOOLS::Blob * NextScatter(const double & bfac=1.);

    inline const double & muF2() const { return m_muf2; }
    inline const double & muR2() const { return m_mur2; }
    
    inline void SetMIProcesses(MI_Processes * processes) {
      p_processes = processes;
    }
    inline void SetOverEstimator(Over_Estimator * overestimator) {
      p_overestimator = overestimator;
    }
    inline void SetResidualX(const double & x1=1., const double & x2=1.) {
      m_residualx1 = x1;
      m_residualx2 = x2;
    }
    inline void SetLastPT2(const double pt2=-1.) {
      m_lastpt2 = pt2<0.?m_S/4.:pt2; 
    }
    inline const double & LastPT2() const { return m_lastpt2; }
    inline const double & PT2()     const { return m_pt2; }
    inline const double & Y3()      const { return m_y3; }
    inline const double & Y4()      const { return m_y4; }

    
    void Test(const double & Q2,const long int & n);
  };
}

#endif
