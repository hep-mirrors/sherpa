#ifndef DIM__Main__MCatNLO_H
#define DIM__Main__MCatNLO_H

#include "PDF/Main/NLOMC_Base.H"

#include "DIM/Main/Gamma.H"
#include "DIM/Shower/Shower.H"
#include "DIM/Tools/Amplitude.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace DIM {

  class MCatNLO: public PDF::NLOMC_Base {
  private:

    Shower *p_mcatnlo;
    Gamma  *p_gamma;

    Amplitude_Vector m_ampls;

    ATOOLS::Mass_Selector     *p_ms;
    ATOOLS::Cluster_Amplitude *p_rampl;

    int m_psmode, m_wcheck;
    double m_maxweight;

    Amplitude *Convert(ATOOLS::Cluster_Amplitude *const campl,
		       std::map<ATOOLS::Cluster_Leg*,Parton*> &lmap);
    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl,
		       const bool & soft=false);
    void CleanUp();

  public:

    MCatNLO(const PDF::NLOMC_Key &key);

    ~MCatNLO();

    int GeneratePoint(ATOOLS::Cluster_Amplitude *const ampl);

    double KT2(const ATOOLS::NLO_subevt &sub,
	       const double &x,const double &y,const double &Q2,
           const double paipb=0, const double alpha=0, const double cos_phi=0);

    ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude(const int mode=0);
    inline ATOOLS::Cluster_Amplitude *GetBornAmplitude() { return p_rampl; }

    inline int PSMode() const { return m_psmode; }

  };// end of class MCatNLO

}// end of namespace DIM

#endif
