#ifndef AMISIC_Tools_Impact_Parameter_H
#define AMISIC_Tools_Impact_Parameter_H

#include "AMISIC++/Tools/Interaction_Probability.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace AMISIC {
  class MI_Processes;
  
  class Impact_Parameter {
    Interaction_Probability * p_pint;
    Matter_Overlap          * p_mo;
    MI_Processes            * p_procs;
    double m_fc, m_oexp, m_enhancement;
    double m_b, m_bmax, m_bnorm;
    bool   m_test, m_ana;
    std::map<std::string,ATOOLS::Histogram *> m_histos;

    void InitAnalysis();
    void BAnalyse(const double & pt2,const double & b);
    void Analyse(const double & pt2,const double & sudakov,
		 const double & softpart, const double & hardpart);
    void FinishAnalysis();
    void Test();
  public:
    Impact_Parameter();
    ~Impact_Parameter();

    void   Initialize(const double & xsecratio);
    void   Update(double xsecratio);
    double SelectB(const double & pt2);
    double operator()(const double & b);

    const double & Enhancement() const { return m_enhancement; }

    inline void SetProcesses(MI_Processes * procs) { p_procs = procs; }
  };
}

#endif
