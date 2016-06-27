#ifndef PHASIC_Process_MCATNLO_Process_H
#define PHASIC_Process_MCATNLO_Process_H

#include "ATOOLS/Org/Smart_Pointer.H"
#include "SHERPA/Tools/Variations.H"

#include "PHASIC++/Process/Process_Base.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace PDF {
  class Cluster_Definitions_Base;
  class NLOMC_Base;
}

namespace MODEL { class Running_AlphaS; }

namespace PHASIC {

  class ME_Generators;
  class MCatNLO_Process;

  class MCatNLO_Process : public Process_Base {
  private:

    ME_Generators &m_gens;

    Process_Base     *p_bviproc, *p_rsproc, *p_bproc, *p_rproc, *p_ddproc;
    PDF::NLOMC_Base  *p_nlomc;

    double m_lastrs, m_lastb, m_lastvi;
    int    m_fomode, m_hpsmode, m_kfacmode, m_wassevent;

    std::string m_rsscale;

    ATOOLS::StringIDipSet_Map m_iinfo;
    ATOOLS::StringDDipSet_Map m_dinfo;

    ATOOLS::Cluster_Amplitude *p_ampl;

    Process_Base *InitProcess(const Process_Info &pi,
			      const nlo_type::code nlotype,const bool real);

    Process_Base *FindProcess(const ATOOLS::Cluster_Amplitude *ampl,
			      const nlo_type::code type=nlo_type::lo,
			      const bool error=true) const;
    Process_Base *FindProcess(const ATOOLS::NLO_subevt *sub,
			      const nlo_type::code type=nlo_type::lo) const;

    ATOOLS::Cluster_Amplitude *CreateAmplitude
    (const ATOOLS::NLO_subevt *sub) const;

    double OneSEvent(const int wmode);
    double OneHEvent(const int wmode);

    /*!
     * @param bvi  double Differential XS for BVI
     * @param b    double Differential XS for B
     * @param rs   double Differential XS for RS
     * @param r    double Differential XS for R 
     * @param ran  double Random number used to select S/H
     * @param ampl Cluster_Amplitude Amplitude which should be updated with the
     *             local K factor when S is selected
     * @return
     *
     * This method calculates the local Kfactor for MENLOPS for given
     * differential XS values
     */
    double LocalKFactor(double bvi, double b, double rs, double r,
                        double random,
                        const ATOOLS::Cluster_Amplitude *ampl=NULL);

    struct KFactorReweightingInfo {
      double m_random;
      SP(SHERPA::Variation_Weights) m_rsvarweights;
      SP(SHERPA::Variation_Weights) m_bvivarweights;
      SP(SHERPA::Variation_Weights) m_bvarweights;
    };

    double ReweightLocalKFactor(SHERPA::Variation_Parameters *,
                                SHERPA::Variation_Weights *,
                                KFactorReweightingInfo &);

  public:

    MCatNLO_Process(ME_Generators &gens,NLOTypeStringProcessMap_Map *pmap);

    ~MCatNLO_Process();

    void Init(const Process_Info &pi,
              BEAM::Beam_Spectra_Handler *const beam,
              PDF::ISR_Handler *const isr,const int mode=0);
    bool InitSubtermInfo();

    bool IsGroup() const;
    size_t Size() const;
    Process_Base *operator[](const size_t &i);

    ATOOLS::Weight_Info *OneEvent(const int wmode,const int mode=0);

    double Differential(const ATOOLS::Vec4D_Vector &p);

    bool CalculateTotalXSec(const std::string &resultpath,
                            const bool create=false);
    void SetLookUp(const bool lookup);
    void SetScale(const Scale_Setter_Arguments &scale);
    void SetKFactor(const KFactor_Setter_Arguments &args);
    void SetFixedScale(const std::vector<double> &s);
    void SetSelector(const Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);
    void SetVariationWeights(SHERPA::Variation_Weights *const);
    void SetClusterDefinitions
    (PDF::Cluster_Definitions_Base *const cluster);

    double LocalKFactor(const ATOOLS::Cluster_Amplitude &ampl);

    ATOOLS::Cluster_Amplitude *GetAmplitude();

    void InitPSHandler(const double &maxerror,
		       const std::string eobs,
		       const std::string efunc);

    inline void SetMCatNLO(PDF::NLOMC_Base *const mc) { p_nlomc=mc; }

    bool WasSEvent() const { return m_wassevent; }

  };// end of class MCatNLO_Process

}// end of namespace PHASIC

#endif
