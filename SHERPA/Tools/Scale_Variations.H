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
    inline double Value() const  { return m_val; }

    inline int    Pdf1Id() const  { return m_pdf1id; }
    inline int    Pdf2Id() const  { return m_pdf2id; }
    inline int    PdfId() const
    { if (m_pdf1id==m_pdf2id) return m_pdf1id; else return -1; }

    inline PDF::PDF_Base * PDF1() { return p_pdf1; }
    inline PDF::PDF_Base * PDF2() { return p_pdf2; }

    inline MODEL::One_Running_AlphaS * AlphaS() { return p_as; }

    inline std::string Name() const { return m_name; }

    inline void SetValue(const double &val) { m_val=val; }

  };

  typedef std::map<std::string, Scale_Variation *> NamedScaleVariationMap;

  std::ostream& operator<<(std::ostream &s,const Scale_Variation &sv);
  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap &nsvm);
  std::ostream& operator<<(std::ostream &s,const NamedScaleVariationMap *nsvm);

  class Scale_Variations {
  private:
    bool m_on,m_loadlhapdf;
    NamedScaleVariationMap * p_nsvmap;

    struct Parameters {
      double B,VI,KP,RS,DADS,renterm,facterm1,facterm2;
      double muR2,muF12,muF22;
      double x1,x2,x1p,x2p,x1DADS,x2DADS;
      std::vector<double> wren,wfac;
      double oqcd,oew;
      ATOOLS::Flavour fl1,fl2;
      ATOOLS::mewgttype::code type;
    } m_params;

    void ResetValues();
    void ExtractParameters(const ATOOLS::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);
    bool Calculate(Scale_Variation * sv);
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
