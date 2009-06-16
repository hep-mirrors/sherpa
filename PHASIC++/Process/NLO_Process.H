#ifndef PHASIC_PROCESS_NLO_PROCESS_H
#define PHASIC_PROCESS_NLO_PROCESS_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generators.H"

namespace ATOOLS {
  class Cluster_Amplitude;
  class Cluster_Definitions_Base;
}

namespace AMEGIC {
  class Single_DipoleTerm;
}

namespace PHASIC {
  
  class NLO_Process : public Process_Base
  {
    ME_Generators& m_gens;
    Process_Base *p_bproc, *p_vproc, *p_iproc, *p_rproc, *p_sproc;

    std::vector<std::vector<size_t> >   m_rbmap;
    std::vector<double>                 m_bxsecs;

    Process_Map m_pmap_born;
    Process_Map m_pmap_real;
    
    Process_Base* InitBornLike(nlo_type::code nlotype);
    Process_Base* InitRealLike(nlo_type::code nlotype);
    ATOOLS::Cluster_Amplitude * CreateAmplitude(const Process_Base * proc,
                                                const ATOOLS::Vec4D_Vector &p);
    void SetIds(ATOOLS::Cluster_Amplitude *bam, ATOOLS::Cluster_Amplitude *ram,
                AMEGIC::Single_DipoleTerm *sdt);

    void FillRBMap();
    Process_Base * FindBorn(const ATOOLS::Cluster_Amplitude * born);
    Process_Base * FindReal(const ATOOLS::Cluster_Amplitude * real);

  public:
    NLO_Process(ME_Generators& gens);

    void Init(const Process_Info &pi,
              BEAM::Beam_Spectra_Handler *const beam,
              PDF::ISR_Handler *const isr);

    size_t Size() const;
    Process_Base *operator[](const size_t &i);

    void DeSelect();
    bool SelectOne();

    Weight_Info *OneEvent();
    Weight_Info *WeightedEvent(const int mode=0);

    double Differential(const ATOOLS::Vec4D_Vector &p);
    double Differential2();
    void AddPoint(const double &value);
    double RealEmissionWeight(const ATOOLS::Cluster_Amplitude * realampl);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);
    void BuildCuts(Cut_Data *const cuts);
    void UpdateCuts(const double &sp,const double &y,
		    Cut_Data *const cuts);

    bool FillIntegrator(Phase_Space_Handler *const psh);
    void UpdateIntegrator(Phase_Space_Handler *const psh);

    bool CalculateTotalXSec(const std::string &resultpath,
                            const bool create=false);
    void SetLookUp(const bool lookup);
    void SetScale(const std::string &scale,
		  const std::string &mur2tag,
		  const std::string &muf2tag);
    void SetKFactor(const std::string &kfactor,
		    const size_t &oqcdlo,const size_t &oewlo);
    void SetKFactorOn(const bool on);
    void SetSelector(const Selector_Key &key);

    void SetClusterDefinitions
    (ATOOLS::Cluster_Definitions_Base *const cluster);

  };

}

#endif // NLO_PROCESS_H
