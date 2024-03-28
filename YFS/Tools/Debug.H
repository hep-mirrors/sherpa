#ifndef YFS_Main_DEBUG_H
#define YFS_Main_DEBUG_H

#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/MathTools.H"

#include "YFS/Main/YFS_Base.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/FSR.H"


#include <map>

namespace YFS {
typedef std::map<std::string, ATOOLS::Histogram *> histograms;

  class Debug : public YFS_Base
  {
  public:
    Debug();
    ~Debug();

    histograms m_histograms_ISR, m_histograms_FSR;

    void InitializeHist();
    void FillHist(const Vec4D_Vector &plab, YFS::ISR *p_isr, YFS::FSR *p_fsr, double weight=1);
    void WriteHistograms();

  };

}
#endif