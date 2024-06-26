#ifndef PHASIC_Process_MCATNLO_Process_H
#define PHASIC_Process_MCATNLO_Process_H

#include "ATOOLS/Phys/Variations.H"

#include "PHASIC++/Process/Process_Base.H"

#include <memory>

namespace ATOOLS { class Cluster_Amplitude; }

namespace PDF {
  class Cluster_Definitions_Base;
}

namespace MODEL { class Running_AlphaS; }

namespace PHASIC {

  class ME_Generators;
  class Jet_Finder;

  //! This class manages all processes contributing to an MC@@NLO calculation
  class MCatNLO_Process : public Process_Base {
  private:

    //! Reference to the matrix element generators
    ME_Generators &m_gens;

    //! Pointers to the contributing processes
    Process_Base     *p_bviproc, *p_rsproc, *p_bproc, *p_rproc, *p_ddproc;
    //! Pointer to the MC@@NLO shower
    PDF::NLOMC_Base  *p_nlomc;

    //! Book-keeping variables for last RS, B, VI contributions
    double m_lastrs, m_lastb, m_lastvi;
    //! Switches for the behaviour of this class
    int    m_fomode, m_psmode, m_hpsmode, m_kfacmode;
    //! Book-keeping to track whether event was S or H
    bool   m_wassevent;
    //! Separate scale for RS events
    std::string m_rsscale;
    
    //! Map of integrated dipoles
    ATOOLS::StringIDipSet_Map m_iinfo;
    //! Map of differential dipoles
    ATOOLS::StringDDipSet_Map m_dinfo;

    //! Cluster amplitude for interface with MC@@NLO shower
    ATOOLS::Cluster_Amplitude *p_ampl;

    void RegisterDefaults() const;

    /*!
     * @param pi Process_Info Process information object
     * @param nlotype nlo_type::code NLO type of process
     * @param real bool Whether this is a real correction or not
     * @return Pointer to the initialised process
     *
     * This method initialises all needed processes.
     */
    Process_Base *InitProcess(const Process_Info &pi,
			      const ATOOLS::nlo_type::code nlotype,
			      const int real);

    /*!
     * @param ampl Cluster_Amplitude Cluster amplitude for which the process is
     *   needed
     * @param type NLO type of process needed
     * @param error bool Exit when not found or not
     * @return Pointer to the initialised process
     *
     * This method retrieves processes from the global process map.
     */
    Process_Base *FindProcess(const ATOOLS::Cluster_Amplitude *ampl,
                              const ATOOLS::nlo_type::code type
                                    =ATOOLS::nlo_type::lo,
			      const bool error=true) const;

    /*!
     * @param sub NLO_subevt NLO subevent for which the cluster amplitude
     *                       should be constructed
     * @return Pointer to cluster amplitude constructed
     *
     * This method constructs a cluster amplitude for the given NLO subevent
     */
    ATOOLS::Cluster_Amplitude *CreateAmplitude
    (const ATOOLS::NLO_subevt *sub) const;

    /*!
     * @param wmode int Mode (not used)
     * @return additional S-Event weights (i.e. weights from the MC@@NLO
     * shower), or zero
     *
     * This method generates one S-event
     */
    ATOOLS::Weights_Map OneSEvent(const int wmode);
    /*!
     * @param wmode int Mode (not used)
     * @return additional H-Event weight (one or zero)
     *
     * This method generates one H-event
     */
    ATOOLS::Weights_Map OneHEvent(const int wmode);

    struct LocalKFactorInfo {
      double s {0.0};
      double h {0.0};
    };
    /*!
     * @param bvi  double Differential XS for BVI
     * @param b    double Differential XS for B
     * @param rs   double Differential XS for RS
     * @param r    double Differential XS for R
     * @return S and H values corresponding to m_kfacmode
     *
     * This method calculates S and H values needed to calculate the local
     * Kfactor for MENLOPS for given differential XS values
     */
    LocalKFactorInfo
    CalculateLocalKFactorInfo(double bvi, double b, double rs, double r);

    struct HEventVeto_Args {
      Jet_Finder *p_jf;
      double m_jcv;
      HEventVeto_Args(Jet_Finder *const jf,const double &jcv):
	p_jf(jf), m_jcv(jcv) {}
    };

  public:

    /*!
     * @param gens ME_Generators Matrix element generators
     * @param pmap NLOTypeStringProcessMap_Map global process map
     *
     * Constructor
     */
    MCatNLO_Process(ME_Generators &gens,NLOTypeStringProcessMap_Map *pmap);

    /*!
     * Destructor
     */
    ~MCatNLO_Process();

    /*!
     * @param pi Process_Info Process information on the desired process
     * @param beam Beam_Spectra_Handler Beam spectra handler
     * @param isr ISR_Handler ISR Handler
     * @param mode int Mode (not used)
     *
     * This method initialises the MC@@NLO process and all contributions
     */
    void Init(const Process_Info &pi,
              BEAM::Beam_Spectra_Handler *const beam,
              PDF::ISR_Handler *const isr, YFS::YFS_Handler *const yfs,
              const int mode=0) override;
    /*!
     * @return Returns true
     *
     * This method initialises m_iinfo and m_dinfo
     */
    bool InitSubtermInfo();

    /*!
     * @return Returns true
     *
     * Inherited from Process_Base, returns whether this process has subprocesses
     */
    bool IsGroup() const override;
    /*!
     * @return Returns 2
     *
     * Inherited from Process_Base, returns the number of subprocesses
     */
    size_t Size() const override;
    /*!
     * @param i size_t Index
     * @return Returns true
     *
     * Inherited from Process_Base, returns p_bviproc (i=0) or p_rsproc (i=1)
     */
    Process_Base *operator[](const size_t &i) override;

    /*!
     * @param wmode int Passed on to methods called
     * @param varmode Variations_Mode Specifies which weights to calculate
     * @param mode int Passed on to methods called
     * @return Returns weight information object
     *
     * Inherited from Process_Base, generates one event
     */
    ATOOLS::Weight_Info *OneEvent(
        const int wmode,
        ATOOLS::Variations_Mode varmode=ATOOLS::Variations_Mode::all,
        const int mode=0) override;

    /*!
     * @param p Four momenta
     * @param varmode Specifies which weights to calculate (all,
     * nominal-only, ...)
     * @return Returns weight of given momentum configuration
     *
     * Inherited from Process_Base, calculates the weight of given configuration
     */
    ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector &,
                                     ATOOLS::Variations_Mode) override;

    /*!
     * @param resultpath string Path where to store the results
     * @param create bool Whether to create the results file or not
     * @return Returns weight information object
     *
     * Inherited from Process_Base, asks subprocesses to integrate themselves
     */
    bool CalculateTotalXSec(const std::string &resultpath,
                            const bool create=false) override;
    /*!
     * @param lookup bool Whether results can be looked up
     *
     * Inherited from Process_Base, set whether mapped processes lookup their
     * result or calculate anew, sets it for all subprocesses
     */
    void SetLookUp(const bool lookup) override;
    /*!
     * Inherited from Process_Base, calls InitScale for all subprocesses
     */
    bool InitScale() override;
    /*!
     * @param scale Scale_Setter_Arguments Scale setter arguments
     *
     * Inherited from Process_Base, calls SetScale for all subprocesses
     */
    void SetScale(const Scale_Setter_Arguments &scale) override;
    /*!
     * @param args KFactor_Setter_Arguments KFactor setter arguments
     *
     * Inherited from Process_Base, calls SetKFactor for all subprocesses
     */
    void SetKFactor(const KFactor_Setter_Arguments &args) override;
    /*!
     * @param mode Variations_Mode
     *
     * Inherited from Process_Base, calls InitializeTheReweighting for all
     * subprocess
     */
    void InitializeTheReweighting(ATOOLS::Variations_Mode) override;
    /*!
     * @param s vector<double> Fixed scales
     *
     * Inherited from Process_Base, calls SetFixedScale for all subprocesses
     */
    void SetFixedScale(const std::vector<double> &s) override;
    /*!
     * @param key Selector_Key Selector key
     *
     * Inherited from Process_Base, calls SetSelector for all subprocesses
     */
    void SetSelector(const Selector_Key &key) override;
    /*!
     * @param ps Shower_Base Parton shower
     *
     * Inherited from Process_Base, calls SetShower for all subprocesses
     */
    void SetShower(PDF::Shower_Base *const ps) override;
    void SetNLOMC(PDF::NLOMC_Base *const mc) override;

    /*!
     * @param ampl Cluster_Amplitude Cluster amplitude (might be updated)
     * @return Local KFactor
     *
     * This method calculates the local Kfactor for MENLOPS for given config
     */
    ATOOLS::Weights_Map LocalKFactor(ATOOLS::Cluster_Amplitude& ampl);

    /*!
     * @return Cluster amplitude p_ampl
     *
     * This method return p_ampl
     */
    ATOOLS::Cluster_Amplitude *GetAmplitude();

    /*!
     * @param maxerror double Error target
     * @param eobs string Enhance observable
     * @param efunc string Enhance function
     *
     * Inherited from Process_Base, calls InitPSHandler for all subprocesses
     */
    void InitPSHandler(const double &maxerror,
		       const std::string eobs,
		       const std::string efunc) override;

    /*!
     * @param mc NLOMC_Base MC@@NLO shower
     *
     * This method sets the MC@@NLO shower
     */
    inline void SetMCatNLO(PDF::NLOMC_Base *const mc) { p_nlomc=mc; }

    /*!
     * @return Whether event last generated was S-Event or not
     *
     * This method returns whether last generated event was S-Event or not
     */
    bool WasSEvent() const { return m_wassevent; }
    Process_Base* BVIProc() { return p_bviproc; }

  };// end of class MCatNLO_Process

}// end of namespace PHASIC

#endif
