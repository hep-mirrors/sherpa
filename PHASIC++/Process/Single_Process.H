#ifndef PHASIC_Process_Single_Process_H
#define PHASIC_Process_Single_Process_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Variations.H"

namespace PHASIC {

  class Single_Process: public Process_Base {

  protected:

    double m_lastxs, m_lastbxs, m_lastflux;
    bool   m_zero, m_pdfcts, m_dads;

    //! Override in a subclass that calculates KP terms to enable their reweighting
    virtual double KPTerms(int mode, double scalefac2=1.0) { return 0.; };

  public:

    Single_Process();
    virtual ~Single_Process();

    size_t Size() const;
    Process_Base *operator[](const size_t &i);

    ATOOLS::Weight_Info *OneEvent(const int wmode,const int mode=0);

    double KFactor(const int mode=0) const;

    double CollinearCounterTerms
    (const int i,const ATOOLS::Flavour &fl,const ATOOLS::Vec4D &p,
     const double &z,const double &t1,const double &t2,
     const double &muf2fac=1.0,
     const double &mur2fac=1.0,
     MODEL::One_Running_AlphaS * as=NULL) const;

    double Differential(const ATOOLS::Vec4D_Vector &p);

    bool CalculateTotalXSec(const std::string &resultpath,
                            const bool create);
    void SetLookUp(const bool lookup);

    void SetScale(const Scale_Setter_Arguments &args);
    void SetKFactor(const KFactor_Setter_Arguments &args);

    ATOOLS::Cluster_Amplitude *Cluster(const ATOOLS::Vec4D_Vector &p,
                                       const size_t &mode=0);

    virtual double Partonic(const ATOOLS::Vec4D_Vector &p,
                            const int mode) = 0;

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
        bool skipsfirstampl, const double &Q2,
        const double &muf2fac=1.0,
        const double &mur2fac=1.0,
        const double &showermuf2fac=1.0,
        MODEL::One_Running_AlphaS * as=NULL,
        const ATOOLS::Cluster_Sequence_Info * const nominalcsi=NULL);
    void AddISR(ATOOLS::Cluster_Sequence_Info &,
                const ATOOLS::ClusterAmplitude_Vector &,
                bool skipsfirstampl, const double &Q2,
                const double &muf2fac=1.0,
                const double &mur2fac=1.0,
                const double &showermuf2fac=1.0,
                MODEL::One_Running_AlphaS * as=NULL,
                const ATOOLS::Cluster_Sequence_Info * const nominalcsi=NULL);
    void AddBeam(ATOOLS::Cluster_Sequence_Info &, const double &Q2);

    // Reweighting utilities and functions
    struct BornLikeReweightingInfo {
      double m_wgt;
      size_t m_orderqcd;
      int m_fl1, m_fl2;
      double m_x1, m_x2;
      double m_muR2, m_muF12, m_muF22;
      ATOOLS::ClusterAmplitude_Vector m_ampls;
      bool m_skipsfirstampl;
      double m_fallbackresult;
    };
    //! Reweighting events with subevents (i.e. RS)
    ATOOLS::Subevent_Weights_Vector ReweightSubevents(
        ATOOLS::Variation_Parameters *,
        ATOOLS::Variation_Weights *,
        const long &additionaldata);
    //! Reweighting events without subevents (i.e. *except* RS)
    double ReweightWithoutSubevents(ATOOLS::Variation_Parameters *,
                                    ATOOLS::Variation_Weights *,
                                    ATOOLS::ClusterAmplitude_Vector &);
    //! Return (wgt * \alpha_S * csi.pdfwgt)
    double ReweightBornLike(ATOOLS::Variation_Parameters *,
                            BornLikeReweightingInfo &);
    //! Return (pdf1, pdf2)
    std::pair<double, double> GetPairOfPDFValuesOrOne(
        ATOOLS::Variation_Parameters *,
        Single_Process::BornLikeReweightingInfo &) const;
    //! Generate cluster sequence info for variation parameters
    ATOOLS::Cluster_Sequence_Info ClusterSequenceInfo(
        ATOOLS::Variation_Parameters * varparams,
        BornLikeReweightingInfo & info,
        const double &mur2fac,
        const ATOOLS::Cluster_Sequence_Info * const nominalcsi);
    //! Generate KP terms weight for variation parameters
    double KPTerms(ATOOLS::Variation_Parameters *);
    //! Calculate MuR for variation parameters
    double MuR2(ATOOLS::Variation_Parameters *,
                Single_Process::BornLikeReweightingInfo &) const;
    //! Calculate AlphaS ratio (*asnew)(mur2new) / MODEL::as(mur2old)
    double AlphaSRatio(double mur2old, double mur2new,
                       MODEL::One_Running_AlphaS * asnew);
  };// end of class Single_Process

}// end of namespace PHASIC

#endif
