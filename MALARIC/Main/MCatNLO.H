#ifndef MALARIC__Main__MCatNLO_H
#define MALARIC__Main__MCatNLO_H

#include "PDF/Main/NLOMC_Base.H"

#include "MALARIC/Main/Gamma.H"
#include "MALARIC/Shower/Shower.H"
#include "MALARIC/Tools/Amplitude.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace MALARIC {

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
	       const PHASIC::Kin_Args *ka,
	       const PHASIC::Ant_Args *aa);

    ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude(const int mode=0);
    inline ATOOLS::Cluster_Amplitude *GetBornAmplitude() { return p_rampl; }

    inline int PSMode() const { return m_psmode; }

  };// end of class MCatNLO

}// end of namespace MALARIC

#endif
