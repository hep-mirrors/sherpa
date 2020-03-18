#ifndef PHASIC_Main_Phase_Space_Enhance_H
#define PHASIC_Main_Phase_Space_Enhance_H

#include "PHASIC++/Enhance/Enhance_Observable_Base.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <iostream>

namespace PHASIC {
  class Process_Base;
  class Phase_Space_Handler;

  class Phase_Space_Enhance {
  private:
    Enhance_Observable_Base * p_obs, * p_func;
    ATOOLS::Histogram       * p_histo, * p_histo_current;    
    int                       m_xs;
    double                    m_factor;
    ATOOLS::Vec4D           * p_moms;
    ATOOLS::Flavour         * p_flavs;

    size_t m_nflavs;

    void RegisterDefaults();
  public:
    Phase_Space_Enhance();
    ~Phase_Space_Enhance();

    void   Init(Phase_Space_Handler * psh);
    double operator()();
    double Factor(Process_Base *const proc,const double & totalxs);

    void SetObservable(const std::string &enhanceobs,Process_Base * const process);
    void SetFunction(const std::string &enhancefunc,Process_Base * const process);
    
    void AddPoint(const double xs,Process_Base * const process);
    void Optimize();
    void ReadIn(const std::string &path);

    inline void MPISync() {
      if (p_histo_current) p_histo_current->MPISync();
    }
    inline void WriteOut(const std::string &path) {
      if (p_histo) p_histo->Output(path+"/MC_Enhance.histo");
    }
  };
}
#endif
