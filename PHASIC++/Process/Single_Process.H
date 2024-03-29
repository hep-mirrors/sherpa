#ifndef PHASIC_Process_Single_Process_H
#define PHASIC_Process_Single_Process_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Variations.H"

namespace ATOOLS {
  class Hard_Process_Variation_Generator_Base;
}

namespace PHASIC {

  class Single_Process: public Process_Base {

  protected:

    double m_lastxs, m_lastbxs, m_dsweight, m_lastflux;
    ATOOLS::Weights_Map m_dadswgtmap;
    bool   m_zero, m_pdfcts, m_dads;
    size_t m_nfconvscheme;
    ATOOLS::Cluster_Sequence_Info m_csi;
    std::vector<ATOOLS::asscontrib::type> m_asscontrib;

    //! Override in a subclass that calculates KP terms to enable their reweighting
    virtual double KPTerms(int mode, double scalefac2=1.0) { return 0.; };

    //! Override in a subclass as a fallback if the underlying ME returns the
    //! special QCD order of 98, thus indicating that it would like to report
    //! its own renormalization scale reweighting factor.
    virtual double CustomRelativeVariationWeightForRenormalizationScaleFactor(double) const {
      THROW(not_implemented, "This must be implemented in a subclass.");
    }

  public:

    Single_Process();
    virtual ~Single_Process();

    size_t Size() const override;
    Process_Base *operator[](const size_t &i) override;

    ATOOLS::Weight_Info *OneEvent(
        const int wmode,
        ATOOLS::Variations_Mode varmode=ATOOLS::Variations_Mode::all,
        const int mode=0) override;

    double KFactor(const int mode=0) const;

    double NfSchemeConversionTerms() const;
    double CollinearCounterTerms
    (const int i,const ATOOLS::Flavour &fl,const ATOOLS::Vec4D &p,
     const double &z,const double &t1,const double &t2,
     const double &muf2fac=1.0,
     const double &mur2fac=1.0,
     MODEL::Running_AlphaS * as=NULL) const;

#if 0
    double Differential(const ATOOLS::Vec4D_Vector &p);
#endif

    bool CalculateTotalXSec(const std::string &resultpath,
                            const bool create) override;
    void SetLookUp(const bool lookup) override;

    void SetScale(const Scale_Setter_Arguments &args) override;
    void SetKFactor(const KFactor_Setter_Arguments &args) override;
    void InitializeTheReweighting(ATOOLS::Variations_Mode) override;

    ATOOLS::Cluster_Amplitude *Cluster(const ATOOLS::Vec4D_Vector &p,
                                       const size_t &mode=0);

    /*!
     * Calculate the partonic cross section
     *
     * This updates m_last, m_lastb, m_lastxs and corresponding fields of
     * m_mewgtinfo. These information can be retrieved using Last, LastB,
     * LastXS and GetMEwgtinfo.
    */
    virtual double Partonic(
        const ATOOLS::Vec4D_Vector &p,
        ATOOLS::Variations_Mode varmode=ATOOLS::Variations_Mode::all,
        int mode=0) = 0;

    inline virtual const double SPrimeMin() const override { return -1.; }
    inline virtual const double SPrimeMax() const override { return -1.; }
    inline virtual const bool   HasInternalScale() const override { return false; }
    inline virtual const double InternalScale()    const override { return -1.; }

    virtual bool Combinable(const size_t &idi,
                            const size_t &idj);
    virtual const ATOOLS::Flavour_Vector &
    CombinedFlavour(const size_t &idij);

    virtual ATOOLS::Flavour ReMap
    (const ATOOLS::Flavour &fl,const size_t &id) const;

    inline void ResetLastXS() { m_lastxs=0.0; }

    inline double LastXS() const { return m_lastxs; }

    inline bool Zero() const { return m_zero; }

  private:

    //! Generate cluster sequence info (containing the ISR+beam weight)
    ATOOLS::Cluster_Sequence_Info ClusterSequenceInfo(
        const ATOOLS::ClusterAmplitude_Vector &,
        const double &Q2,
        const double &muf2fac=1.0,
        const double &mur2fac=1.0,
        const double &showermuf2fac=1.0,
        MODEL::Running_AlphaS * as=NULL,
        const ATOOLS::Cluster_Sequence_Info * const nominalcsi=NULL);
    void AddISR(ATOOLS::Cluster_Sequence_Info &,
                const ATOOLS::ClusterAmplitude_Vector &,
                const double &Q2,
                const double &muf2fac=1.0,
                const double &mur2fac=1.0,
                const double &showermuf2fac=1.0,
                MODEL::Running_AlphaS * as=NULL,
                const ATOOLS::Cluster_Sequence_Info * const nominalcsi=NULL);
    void AddBeam(ATOOLS::Cluster_Sequence_Info &, const double &Q2);

    // Reweighting utilities and functions
    struct BornLikeReweightingInfo {
      BornLikeReweightingInfo() {}
      BornLikeReweightingInfo(const ATOOLS::ME_Weight_Info& mewgtinfo,
                              const ATOOLS::ClusterAmplitude_Vector& ampls,
                              double fallbackresult):
        m_wgt {mewgtinfo.m_B + mewgtinfo.m_VI + mewgtinfo.m_KP},
        m_orderqcd {mewgtinfo.m_oqcd},
        m_fl1 {mewgtinfo.m_fl1},
        m_fl2 {mewgtinfo.m_fl2},
        m_x1 {mewgtinfo.m_x1},
        m_x2 {mewgtinfo.m_x2},
        m_muF2 {mewgtinfo.m_muf2},
        m_muR2 {mewgtinfo.m_mur2},
        m_ampls {ampls},
        m_fallbackresult {fallbackresult}
      {}
      double m_wgt;
      size_t m_orderqcd;
      int m_fl1, m_fl2;
      double m_x1, m_x2;
      double m_muR2, m_muF2;
      ATOOLS::ClusterAmplitude_Vector m_ampls;
      double m_fallbackresult;
    };
    //! Return \f$wgt * \alpha_{S} * csi.pdfwgt\f$
    double ReweightBornLike(ATOOLS::QCD_Variation_Params&,
                            BornLikeReweightingInfo&);
    //! Generate cluster sequence info for variation parameters
    ATOOLS::Cluster_Sequence_Info
    ClusterSequenceInfo(ATOOLS::QCD_Variation_Params& varparams,
                        BornLikeReweightingInfo& info,
                        const double& mur2fac,
                        const ATOOLS::Cluster_Sequence_Info* const nominalcsi);
    //! Generate KP terms weight for variation parameters
    double KPTerms(const ATOOLS::QCD_Variation_Params *);
    virtual double KPTerms(int mode, PDF::PDF_Base *pdfa,
			   PDF::PDF_Base *pdfb, double scalefac2=0.);
    //! Calculate MuR for variation parameters
    double MuR2(const ATOOLS::QCD_Variation_Params&,
                Single_Process::BornLikeReweightingInfo&) const;
    //! Calculate AlphaS ratio (*asnew)(mur2new) / MODEL::as(mur2old)
    double AlphaSRatio(double mur2old, double mur2new,
                       MODEL::Running_AlphaS * asnew);

    ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector&,
                                     ATOOLS::Variations_Mode) override;

    // the following are helper methods for the implementation of Differential
    void ResetResultsForDifferential(ATOOLS::Variations_Mode);
    void InitMEWeightInfo();
    void UpdateIntegratorMomenta(const ATOOLS::Vec4D_Vector&);
    void UpdateIntegratorMomenta(ATOOLS::ClusterAmplitude_Vector&);
    void UpdateSubeventMomenta(ATOOLS::NLO_subevt&);
    void CalculateFlux(const ATOOLS::Vec4D_Vector&);
    void UpdateMEWeightInfo(Scale_Setter_Base*);
    void ReweightBVI(ATOOLS::ClusterAmplitude_Vector&);
    void ReweightRS(ATOOLS::ClusterAmplitude_Vector&);
    void CalculateAssociatedContributionVariations();

    std::vector<ATOOLS::Hard_Process_Variation_Generator_Base*>
      m_hard_process_variation_generators;

  };// end of class Single_Process

}// end of namespace PHASIC

#endif
