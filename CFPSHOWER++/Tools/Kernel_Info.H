#ifndef CFPSHOWER_Tools_Kernel_Info_H
#define CFPSHOWER_Tools_Kernel_Info_H

#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Flavour.H"


namespace CFPSHOWER {
  struct log_type {
    enum code {
      soft = 1,
      coll = 2,
      none = 0
    };
  };
  std::ostream &operator<<(std::ostream &s,const log_type::code & type);

  struct kernel_type {
    enum code {
      FF   = 1,
      FI   = 2,
      IF   = 3,
      II   = 4,
      none = 0
    };
  };
  std::ostream &operator<<(std::ostream &s,const kernel_type::code & type);
  inline kernel_type::code GetCode(const bool & beam1,const bool & beam2) {
    if (beam1>0) return (beam2>0)?kernel_type::II:kernel_type::IF;
    return (beam2>0)?kernel_type::FI:kernel_type::FF;
  };
  
  struct kin_type {
    enum code {
      CataniSeymour = 2, 
      Alaric        = 1,
      none          = 0
    };
  };
  std::ostream &operator<<(std::ostream &s,const kin_type::code & type);

  struct muR_scheme {
    enum code {
      KT2_pipj_for_gqq = 11,
      KT2_all          = 10,
      T_pipj_for_gqq   = 1,
      T_all            = 0
    };
  };
  std::ostream &operator<<(std::ostream &s,const muR_scheme::code & type);

  class Kernel_Info {
    ATOOLS::Flavour           m_split;
    ATOOLS::Flavour_Vector    m_flavs;
    std::vector<size_t>       m_tagsequence;
    log_type::code            m_logtype;
    kin_type::code            m_kintype;
    muR_scheme::code          m_muRscheme;
    kernel_type::code         m_type;
    MODEL::Running_AlphaS   * p_alphaS;
    MODEL::Running_AlphaQED * p_alpha;
    int    m_SFtype, m_kfactor, m_softcorr, m_endpoint, m_cplscheme;            
    double m_asfactor, m_muR2factor;
    void AdjustFlavours();
    void AdjustFlavours12();
    void AdjustFlavours13(ATOOLS::Flavour_Vector & flavs);
  public:
    Kernel_Info() {}
    Kernel_Info(ATOOLS::Flavour & split,
		ATOOLS::Flavour_Vector & flavs,
		const std::vector<size_t> & tagsequence,
	        log_type::code logtype,
		kin_type::code kintype,
		muR_scheme::code muRscheme,
		kernel_type::code type);
    Kernel_Info(log_type::code logtype,
		kin_type::code kintype,
		kernel_type::code type);

    
    inline void SetAlphaS(MODEL::Running_AlphaS * alphaS)    { p_alphaS     = alphaS; }
    inline void SetAlphaQED(MODEL::Running_AlphaQED * alpha) { p_alpha      = alpha; }
    inline void SetKFactor(const int & kfactor)              { m_kfactor    = kfactor; }
    inline void SetSoftCorrection(const int & softcorr)      { m_softcorr   = softcorr; }
    inline void SetSFType(const int & SFtype)                { m_SFtype     = SFtype; }
    inline void SetEndpoint(const int & endpoint)            { m_endpoint   = endpoint; }
    inline void SetCplScheme(const int & cplscheme)          { m_cplscheme  = cplscheme; }
    inline void SetAsFactor(const double & asfactor)         { m_asfactor   = asfactor; }
    inline void SetMuR2Factor(const double & muR2factor)     { m_muR2factor = muR2factor; }

    inline const ATOOLS::Flavour          GetSplit()       const { return m_split; }
    inline const ATOOLS::Flavour_Vector & GetFlavs()       const { return m_flavs; }
    inline const log_type::code           LogType()        const { return m_logtype; } 
    inline const kin_type::code           KinType()        const { return m_kintype; } 
    inline const muR_scheme::code         MuRScheme()      const { return m_muRscheme; } 
    inline const kernel_type::code        Type()           const { return m_type; } 
    inline const std::vector<size_t>    & TagSequence()    const { return m_tagsequence; }
    inline MODEL::Running_AlphaS        * GetAlphaS()      const { return p_alphaS; }
    inline MODEL::Running_AlphaQED      * GetAlpha()       const { return p_alpha; }
    inline const int                    & SFType()         const { return m_SFtype; } 
    inline const int                    & KFactor()        const { return m_kfactor; }
    inline const int                    & Endpoint()       const { return m_endpoint; }
    inline const int                    & SoftCorrection() const { return m_softcorr; }
    inline const int                    & CplScheme()      const { return m_cplscheme; }
    inline const double                 & AsFactor()       const { return m_asfactor; }
    inline const double                 & MuR2Factor()     const { return m_muR2factor; }

    const std::string SFName()  const;
    const std::string GPName()  const;
    const std::string KinName() const;
    const std::string Name()    const;

    friend std::ostream &operator<<(std::ostream &s,const Kernel_Info & info);
  };// end of struct Splitting

  std::ostream &operator<<(std::ostream &s,const Kernel_Info & info);

  inline bool operator<(std::multiset<ATOOLS::Flavour> & ms1,
			std::multiset<ATOOLS::Flavour> & ms2) {
    size_t ms1size(ms1.size()), ms2size(ms2.size());
    if (ms1size<ms2size) return true;
    if (ms1size>ms2size) return false;
    std::multiset<ATOOLS::Flavour>::iterator ms1it=ms1.begin(), ms2it=ms2.begin();
    while (ms1it!=ms1.end() && ms2it!=ms2.end()) {
      if ((*ms1it)<(*ms2it)) return true;
      if ((*ms1it)>(*ms2it)) return false;
      ms1it++;ms2it++;
    }
  };
}// end of namespace CFPSHOWER

#endif
