#ifndef AMISIC_Tools_Impact_Parameter_H
#define AMISIC_Tools_Impact_Parameter_H

#include "AMISIC++/Tools/Interaction_Probability.H"

namespace AMISIC {
  class MI_Processes;
  
  class Impact_Parameter {
    Interaction_Probability * p_pint;
    Matter_Overlap          * p_mo;
    MI_Processes            * p_procs;
    double m_fc, m_oexp, m_enhancement;
    double m_b, m_bmax;
    bool   m_test;
    
    void Test();
  public:
    Impact_Parameter();
    ~Impact_Parameter() { delete p_pint; }

    void   Initialize(const double & xsecratio);
    double SelectB(const double & pt2);
    double operator()(const double & b);

    const double & Enhancement() const { return m_enhancement; }

    inline void SetProcesses(MI_Processes * procs) { p_procs = procs; }
  };
}

#endif
