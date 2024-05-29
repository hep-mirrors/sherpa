#ifndef AMISIC_Tools_Impact_Parameter_H
#define AMISIC_Tools_Impact_Parameter_H

#include "AMISIC++/Tools/Interaction_Probability.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>
#include <string>

namespace AMISIC {
  class MI_Processes;
  class MI_Integrator;
  class Sudakov_Argument;
  
  class Impact_Parameter {
    Interaction_Probability m_pint;
    Matter_Overlap        * p_mo;
    Sudakov_Argument      * p_sudarg;
    MI_Integrator         * p_integrator;

    double FirstB(const double & s);
    double CalculateB(const double & s,const double & scale=-1.);

    bool   m_test, m_ana;
    std::map<std::string,ATOOLS::Histogram *> m_histos;
    void InitAnalysis();
    void BAnalyse(const double & pt2,const double & b);
    void Analyse(const double & pt2,const double & sudarg,
		 const double & softpart, const double & hardpart);
    void FinishAnalysis();
    void Test();
  public:
    Impact_Parameter();
    ~Impact_Parameter();

    void   Initialize(Matter_Overlap * mo,Sudakov_Argument * sud,axis * sbins);
    double operator()(const double & s,const double & scale=-1.);
    double BEnhancement(const double & s,const double & b);

    //////////////////////////////////////////////////////////////////////////
    // Selecting a position for the parton-parton scattering to take place in
    // impact parameter space, once the overall impact parameter has been
    // chosen.
    //////////////////////////////////////////////////////////////////////////
    inline ATOOLS::Vec4D  SelectPositionForScatter(const double & b) const {
      return p_mo->SelectPositionForScatter(b);
    }
    inline const double fb(const double & s,const double & b) {
      return m_pint.fb(s,b);
    }
  };
}

#endif
