#ifndef AMISIC_Tools_Sudakov_Argument_H
#define AMISIC_Tools_Sudakov_Argument_H

#include "AMISIC++/Tools/Lookup_Tables.H"

namespace AMISIC {
  class MI_Processes;

  class Sudakov_Argument {
  private:
    bool             m_variable_s, m_variable_b;
    axis           * p_sbins, * p_bbins, * p_pt2bins;
    OneDim_Table   * p_integral1, * p_function1;
    TwoDim_Table   * p_integral2, * p_function2;
    ThreeDim_Table * p_integral3, * p_function3;
    MI_Processes   * p_processes;
    double           m_norm, m_bJac;
    size_t           m_sbin, m_bbin, m_n_sbins, m_n_bbins;
    bool             m_withOverlap, m_test;

    void   InitTables();
    void   FillTables();
    double CalculateSigmaHard();
    double FillPT2Values(const double & s);
    void   FillBins(const size_t & pt2bin,const double & sigma,
		    const double & integral);
    void   OutputBin(double & integral, const size_t & pt2bin=0);
    void   OutputTables();
  public:
    Sudakov_Argument(MI_Processes * procs);
    ~Sudakov_Argument();

    void   Initialize(axis * sbins,axis * bbins);
    double TrialEvent(const double & s);

    inline axis         * GetPT2Bins()   const { return p_pt2bins; } 
    inline MI_Processes * GetProcesses() const { return p_processes; }
    inline const double operator()(const double & pt2) const
    { return (*p_function1)(pt2); }
    inline const double operator()(const double & s_or_b,
				   const double & pt2) const
    { return (*p_function2)(s_or_b,pt2); }
    inline const double operator()(const double & s,const double & b,
			     const double & pt2) const
    { return (*p_function3)(s,b,pt2); }
  };
}

#endif
