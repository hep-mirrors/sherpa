#ifndef AMISIC_Tools_Sudakov_Argument_H
#define AMISIC_Tools_Sudakov_Argument_H

#include "AMISIC++/Tools/Lookup_Tables.H"

namespace AMISIC {
  class MI_Processes;

  class Sudakov_Argument {
  private:
    axis           m_sbins, m_pt2bins;
    TwoDim_Table   m_integral, m_function;
    MI_Processes * p_processes;
    bool           m_test;

    void FillTables();
    void FillPT2Values(const size_t & sbin,const double & norm);
    void OutputTables();
  public:
    Sudakov_Argument(MI_Processes * procs,const axis & sbins, const axis & pt2bins);

    double XSratio(const double & s);

    inline const axis & GetSbins() const   { return m_sbins; }
    inline const axis & GetPT2bins() const { return m_pt2bins; }
    inline double operator()(const double & s,const double & pt2) const {
      return m_function(s,pt2);
    }
  };
}

#endif
