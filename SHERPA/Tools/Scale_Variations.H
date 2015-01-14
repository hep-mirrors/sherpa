#ifndef SHERPA_Tools_Scale_Variations_H
#define SHERPA_Tools_Scale_Variations_H

#include <map>
#include <string>
#include <vector>
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {
  class ME_wgtinfo;
}

namespace PDF {
  class PDF_Base;
}

namespace PHASIC {
  struct Weight_Info;
  class Process_Base;
}

namespace SHERPA {

  class Scale_Variation {
  private:
    bool m_deletepdfs;
    double m_muR2fac,m_muF2fac,m_val;
    PDF::PDF_Base *p_pdf1,*p_pdf2;
    int m_pdf1id,m_pdf2id;
    std::string m_pdf1set,m_pdf2set;
    int m_pdf1setmember,m_pdf2setmember;
    std::string m_name;

    std::string GenerateName();

  public:
    Scale_Variation(const double &muR2fac, const double &muF2fac,
                    PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                    bool deletepdfs);
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
      double B,V,I,DADS,renterm,facterm1,facterm2;
      double muR2,muF12,muF22;
      double x1,x2,x1p,x2p,x1DADS,x2DADS;
      double oqcd,oew;
      std::vector<double> kpwgts;
      ATOOLS::Flavour fl1,fl2;
    } m_params;

    void ResetValues();
    void ExtractParameters(const PHASIC::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);
    bool Calculate(Scale_Variation * sv);
    double Calculate(const double& B, const double& V, const double& I,
                     const std::vector<double>& wgts,
                     const ATOOLS::Flavour& fl1, const ATOOLS::Flavour& fl2,
                     PDF::PDF_Base * pdf1, PDF::PDF_Base * pdf2,
                     const double& x1, const double& x2,
                     const double& x1p, const double& x2p,
                     const double& muR2,
                     const double& muF12, const double& muF22,
                     const double& muR2fac,const double& muF2fac,
                     const double& oqcd, const double& oew,
                     const size_t& mode);
  public:
    Scale_Variations();
    ~Scale_Variations();

    bool ComputeVariations(const PHASIC::Weight_Info &winfo,
                           PHASIC::Process_Base * proc);

    inline NamedScaleVariationMap * GetNamedScalesMap() const
    { return p_nsvmap; }

  };

  std::ostream& operator<<(std::ostream &s,const Scale_Variations &svs);
}
#endif
