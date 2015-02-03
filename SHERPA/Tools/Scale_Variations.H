#ifndef SHERPA_Tools_Scale_Variations_H
#define SHERPA_Tools_Scale_Variations_H

#include <map>
#include <string>
#include <vector>
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Weight_Info.H"

namespace MODEL {
  class One_Running_AlphaS;
}

namespace PDF {
  class PDF_Base;
}

namespace PHASIC {
  class Process_Base;
}

namespace SHERPA {

  class Scale_Variation {
  private:
    bool m_deletepdfs,m_deleteas;
    double m_muR2fac,m_muF2fac,m_val;
    std::vector<double> m_RSvals;
    PDF::PDF_Base *p_pdf1,*p_pdf2;
    int m_pdf1id,m_pdf2id;
    std::string m_pdf1set,m_pdf2set;
    int m_pdf1setmember,m_pdf2setmember;
    MODEL::One_Running_AlphaS *p_as;
    std::string m_name;

    std::string GenerateName();

  public:
    Scale_Variation(const double &muR2fac, const double &muF2fac,
                    PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                    MODEL::One_Running_AlphaS * as,
                    bool deletepdfs, bool deleteas);
    ~Scale_Variation();

    inline double MuR2Fac() const { return m_muR2fac; }
    inline double MuF2Fac() const { return m_muF2fac; }

    inline double Value() const                 { return m_val; }
    inline double Value(const size_t& i) const  { return m_RSvals[i]; }

    inline int    Pdf1Id() const  { return m_pdf1id; }
    inline int    Pdf2Id() const  { return m_pdf2id; }
    inline int    PdfId() const
    { if (m_pdf1id==m_pdf2id) return m_pdf1id; else return -1; }

    inline PDF::PDF_Base * PDF1() { return p_pdf1; }
    inline PDF::PDF_Base * PDF2() { return p_pdf2; }

    inline MODEL::One_Running_AlphaS * AlphaS() { return p_as; }

    inline std::string Name() const { return m_name; }

    inline void SetValue(const double &val)                { m_val=val; }
    inline void SetValue(const size_t& i,const double &vl) { m_RSvals[i]=vl; }

    inline void DeleteRSValues()                   { m_RSvals.clear(); }
    inline void InitialisRSValues(const size_t& n) { m_RSvals.resize(n,0.); }
  };

  typedef std::map<std::string, Scale_Variation *> NamedScaleVariationMap;

  std::ostream& operator<<(std::ostream &s,const Scale_Variation &sv);
  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap &nsvm);
  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap *nsvm);

  class Scale_Variations {
  private:
    bool m_on,m_loadlhapdf;
    ATOOLS::Flavour m_gluon,m_quark;
    NamedScaleVariationMap * p_nsvmap;

    struct Parameters {
      double B,VI,KP,renterm,facterm1,facterm2;
      double muR2,muF12,muF22;
      double x1,x2,x1p,x2p;
      std::vector<double> rswgts,rsmuR2s,rsmuF2s;
      std::vector<double> renwgts,kpwgts;
      size_t oqcd,oew;
      int fl1,fl2;
      ATOOLS::mewgttype::code type;
      struct DADS {
        double wgt;
        double muR2,muF12,muF22;
        double x1,x2;
        int fl1,fl2;
      };
      std::vector<DADS> dads;
    } m_params;

    void ResetValues();
    void ExtractParameters(const ATOOLS::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);
    bool Calculate(Scale_Variation * sv,
                   PHASIC::Process_Base * proc);
    double Calculate(const double& B, const double& VI,
                     const double& x1, const double& x2,
                     const double& x1p, const double& x2p,
                     const double& muR2,
                     const double& muF12, const double& muF22,
                     const double& muR2fac, const double& muF2fac,
                     const size_t& oqcd, const size_t& oew,
                     const int& fl1, const int& fl2,
                     PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                     MODEL::One_Running_AlphaS * as,
                     const std::vector<double>& renwgts,
                     const std::vector<double>& kpwgts,
                     const ATOOLS::mewgttype::code type);
  public:
    Scale_Variations();
    ~Scale_Variations();

    bool ComputeVariations(const ATOOLS::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);

    inline NamedScaleVariationMap * GetNamedScalesMap() const
    { return p_nsvmap; }

  };

  std::ostream& operator<<(std::ostream &s,const Scale_Variations &svs);
}
#endif
