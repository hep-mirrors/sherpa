#ifndef PHASIC_Process_YFS_Process_H
#define PHASIC_Process_YFS_Process_H



#include "PHASIC++/Process/Process_Base.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "ATOOLS/Phys/Variations.H"
namespace ATOOLS {
  class Hard_Process_Variation_Generator_Base;
}


namespace MODEL {
  class Single_Vertex;
  typedef std::vector<Single_Vertex *> Vertex_List;
}


namespace PHASIC {

  class ME_Generators;

  class YFS_Process : public Process_Base {
  private:

    //! Reference to the matrix element generators
    ME_Generators &m_gens;


    //! Book-keeping variables for last RS, B, VI contributions
    double m_lastrs, m_lastb, m_lastvi;
    //! Switches for the behaviour of this class
    int    m_fomode, m_hpsmode, m_kfacmode;
    //! Book-keeping to track whether event was S or H
    bool   m_wassevent;
    //! Separate scale for RS events
    std::string m_rsscale;
    ME_Generator_Base *p_gen;
    
    //! Map of integrated dipoles
    ATOOLS::StringIDipSet_Map m_iinfo;
    //! Map of differential dipoles
    ATOOLS::StringDDipSet_Map m_dinfo;

    //! Cluster amplitude for interface with YFS
    ATOOLS::Cluster_Amplitude *p_ampl, *p_rampl;
      //! Generate cluster sequence info (containing the ISR+beam weight)
    ATOOLS::Cluster_Sequence_Info ClusterSequenceInfo(
        const ATOOLS::ClusterAmplitude_Vector &,
        const double &Q2,
        const double &muf2fac=1.0,
        const double &mur2fac=1.0,
        const double &showermuf2fac=1.0,
        MODEL::Running_AlphaS * as=NULL,
        const ATOOLS::Cluster_Sequence_Info * const nominalcsi=NULL);
 
    public:
    YFS::YFS_Handler  *p_yfs;
    YFS_Process(ME_Generators &gens,NLOTypeStringProcessMap_Map *pmap);
    void Init(const Process_Info &pi,
              BEAM::Beam_Spectra_Handler *const beam,
              PDF::ISR_Handler *const isr,
              YFS::YFS_Handler *const yfs, const int mode=0) override;

    ~YFS_Process();
    void RegisterDefaults();
    int m_qedmode;
    Process_Base     *p_bornproc, *p_realproc, *p_virtproc, *p_realvirtproc, *p_realrealproc;
    Process_Base* Selected();
    /*!
     * @return Returns true
     *
     * This method initialises m_iinfo and m_dinfo
     */
    bool InitSubtermInfo();
    void FindResonances();
    void FindProcessPossibleResonances(const Flavour_Vector& fv, MODEL::Vertex_List& vlist);
    void InitPSHandler(const double &maxerror,
               const std::string eobs,
               const std::string efunc) override;
    void OneRealEvent();
    inline Process_Base* GetBornProcess() {return p_bornproc;}
    inline Process_Base* GetRealProcess() {return p_realproc;}
    inline Process_Base* GetVirtProcess() {return p_virtproc;}
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
     * @param mode int Passed on to methods called
     * @return Returns weight information object
     *
     * Inherited from Process_Base, generates one event
     */
    ATOOLS::Weight_Info *OneEvent(const int wmode, ATOOLS::Variations_Mode varmode, const int mode=0) override;

    /*!
     * @param p Vec4D_Vector Four momenta
     * @return Returns weight of given momentum configuration
     *
     * Inherited from Process_Base, calculates the weight of given configuration
     // */
    ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector&,
                                     ATOOLS::Variations_Mode) override;
    
    double m_lastxs, m_lastbxs, m_dsweight, m_lastflux, m_lastyfs;

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
    void SetGenerator(ME_Generator_Base *const gen) override;
    inline ME_Generator_Base * Generator() { PRINT_INFO("HERE"); return p_gen;}
    void InitializeTheReweighting(ATOOLS::Variations_Mode mode) override;
    std::vector<ATOOLS::Hard_Process_Variation_Generator_Base*>
      m_hard_process_variation_generators;
    // double Partonic(const ATOOLS::Vec4D_Vector &p, int);
    /*!
     * @param ampl Cluster_Amplitude Cluster amplitude
     * @return Local KFactor
     *
     * This method calculates the local Kfactor for MENLOPS for given config
     */
    double LocalKFactor(const ATOOLS::Cluster_Amplitude &ampl);
    Process_Base* InitProcess(const Process_Info &pi,ATOOLS::nlo_type::code nlotype,const bool real);
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
     * @return Cluster amplitude p_ampl
     *
     * This method return p_ampl
     */
    ATOOLS::Cluster_Amplitude *GetAmplitude();
    };
}

#endif