#ifndef AMISIC_Tools_Impact_Parameter_H
#define AMISIC_Tools_Impact_Parameter_H

#include "AMISIC++/Tools/Interaction_Probability.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace AMISIC {
  class MI_Processes;
  class Sudakov_Argument;
  
  class Impact_Parameter {
    Interaction_Probability * p_pint;
    Matter_Overlap          * p_mo;
    Sudakov_Argument        * p_sudakov;
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

    void   Initialize(MI_Processes * processes);
    void   Update(const double & s);
    double SelectB(const double & s,const double & pt2=-1.);
    double operator()(const double & b);

    inline ATOOLS::Vec4D  SelectPositionForScatter(const double & b) const {
      return p_mo->SelectPositionForScatter(b);
    }
    inline ATOOLS::Vec4D  SelectRelativePositionForParton() const {
      return p_mo->SelectRelativePositionForParton();
    }
  };
}

#endif
