#ifndef ATOOLS_Phys_Weight_Info_H
#define ATOOLS_Phys_Weight_Info_H

#include "ATOOLS/Math/MathTools.H"
#include <ostream>
#include <vector>

namespace ATOOLS {

  struct PDF_Info {
    int m_fl1, m_fl2;
    double m_x1, m_x2, m_muf12, m_muf22, m_xf1, m_xf2;
    PDF_Info(const int &fl1=0, const int &fl2=0,
             const double &x1=0.,const double &x2=0.,
             const double &muf12=0.,const double &muf22=0.,
             const double &xf1=0.,const double &xf2=0.) :
      m_fl1(fl1), m_fl2(fl2), m_x1(x1), m_x2(x2),
      m_muf12(muf12), m_muf22(muf22), m_xf1(xf1), m_xf2(xf2) {}
  }; // end of struct PDF_Info

  inline std::ostream & operator<<(std::ostream & s,
                                   const ATOOLS::PDF_Info & pi)
  {
    s<<" pdf1 = ("<<pi.m_fl1<<","<<pi.m_x1<<","<<sqrt(pi.m_muf12)
                  <<":"<<pi.m_xf1<<") , "
     <<" pdf2 = ("<<pi.m_fl2<<","<<pi.m_x2<<","<<sqrt(pi.m_muf22)
                  <<":"<<pi.m_xf2<<")";
    return s;
  }

  struct Weight_Info {
    PDF_Info m_pdf;
    double m_weight, m_dxs, m_ntrial;
    Weight_Info(const double &weight=0.0,
                const double &dxs=0.0,const double &ntrial=0.0,
                const int &fl1=0, const int &fl2=0,
                const double &x1=0.,const double &x2=0.,
                const double &xf1=0.0,const double &xf2=0.0,
                const double &muf12=0.0,const double &muf22=0.0):
      m_weight(weight), m_dxs(dxs), m_ntrial(ntrial),
      m_pdf(PDF_Info(fl1,fl2,x1,x2,muf12,muf22,xf1,xf2)) {}
  }; // end of struct Weight_Info

  inline std::ostream & operator<<(std::ostream & s,
                                   const ATOOLS::Weight_Info & wi)
  {
    return s<<" w = "<<wi.m_weight<<std::endl;
  }

  struct mewgttype {
    enum code {
      none = 0,
      B    = 1,
      VI   = 2,
      KP   = 4,
      DADS = 8,
      METS = 16,
      H    = 32,
      RS   = 64
    };
  };

  inline mewgttype::code operator|(const mewgttype::code c1,
                                   const mewgttype::code c2)
  { return (mewgttype::code)((int)c1|(int)c2); }
  inline const mewgttype::code &operator|=(mewgttype::code &c1,
                                           const mewgttype::code c2)
  { return c1=(mewgttype::code)((int)c1|(int)c2); }
  inline mewgttype::code operator&(const mewgttype::code c1,
                                   const mewgttype::code c2)
  { return (mewgttype::code)((int)c1&(int)c2); }
  inline const mewgttype::code &operator&=(mewgttype::code &c1,
                                           const mewgttype::code c2)
  { return c1=(mewgttype::code)((int)c1&(int)c2); }
  inline mewgttype::code operator^(const mewgttype::code c1,
                                   const mewgttype::code c2)
  { return (mewgttype::code)((int)c1^(int)c2); }
  inline const mewgttype::code &operator^=(mewgttype::code &c1,
                                           const mewgttype::code c2)
  { return c1=(mewgttype::code)((int)c1^(int)c2); }

  std::ostream & operator<<(std::ostream & s,const mewgttype::code & type);

  struct Cluster_Step {
  public:
    double m_t,m_xa,m_xb,m_xap,m_xbp;
    int m_fla,m_flb;
    double m_pdfrationumerator, m_pdfratiodenominator;
    Cluster_Step(const double& t, const double& xa, const double& xb,
                 const int& fla, const int& flb) :
      m_t(t), m_xa(xa), m_xb(xb), m_xap(1.), m_xbp(1.), m_fla(fla), m_flb(flb),
      m_pdfrationumerator(1.), m_pdfratiodenominator(1.) {}
  };

  inline std::ostream & operator<<(std::ostream & s,
                                   const ATOOLS::Cluster_Step & cs)
  {
    return s<<"("<<cs.m_t<<";"<<cs.m_fla<<","<<cs.m_xa<<"|"
                              <<cs.m_flb<<","<<cs.m_xb<<")";
  }

  struct Cluster_Sequence_Info {
  public:
    double m_pdfwgt, m_flux, m_ct;
    std::vector<Cluster_Step> m_txfl;

    Cluster_Sequence_Info(const double &w=1.0, const double &f=1.0,
                          const double &c=0.0, const double &r=0.0) :
      m_pdfwgt(w), m_flux(f), m_ct(c) {}
    inline void AddSplitting(const double& t,
                             const double& xa, const double& xb,
                             const int& fla, const int& flb)
    { m_txfl.push_back(Cluster_Step(t,xa,xb,fla,flb)); }
    inline void AddFlux(const double& f)        { m_flux*=f; }
    inline void AddWeight(const double& w)      { m_pdfwgt*=w; }
    inline void AddCounterTerm(const double& c,
                               const double& xp, size_t i)
    { m_ct+=c;
      if (i==0) m_txfl.back().m_xap=xp;
      if (i==1) m_txfl.back().m_xbp=xp; }
    inline void AddPDFRatio(const double& n,const double& d)
    { m_txfl.back().m_pdfrationumerator=n;
      m_txfl.back().m_pdfratiodenominator=d; }
  };

  std::ostream & operator<<(std::ostream & s,
                            const ATOOLS::Cluster_Sequence_Info & csi);

  struct DADS_Info {
  public:
    double   m_wgt, m_x1, m_x2;
    int      m_fl1, m_fl2;
    DADS_Info(const double &wgt,
              const double &x1,const double &x2,
              const int &fl1, const int &fl2):
    m_wgt(wgt), m_x1(x1), m_x2(x2), m_fl1(fl1), m_fl2(fl2) {}
  };

  inline std::ostream & operator<<(std::ostream & s,
                                   const ATOOLS::DADS_Info & dadsi)
  {
    s<<"DADS: wgt="<<dadsi.m_wgt<<", x1="<<dadsi.m_x1<<", x2="<<dadsi.m_x2
     <<", fl1="<<dadsi.m_fl1<<", fl2="<<dadsi.m_fl2;
    return s;
  }

  struct RDA_Info {
  public:
    double m_wgt, m_mur2, m_muf12, m_muf22;
    size_t m_i,m_j,m_k;
    RDA_Info(const double &wgt,const double &mur2,
             const double &muf12,const double &muf22,
             const size_t &i,const size_t &j, const size_t &k):
    m_wgt(wgt), m_mur2(mur2), m_muf12(muf12), m_muf22(muf22),
    m_i(i) ,m_j(j) ,m_k(k) {}
  };

  inline std::ostream & operator<<(std::ostream & s,
                                   const ATOOLS::RDA_Info & rdai)
  {
    s<<"RDA("<<rdai.m_i<<","<<rdai.m_j<<","<<rdai.m_k<<"): "
     <<"wgt="<<rdai.m_wgt<<", muR2="<<rdai.m_mur2
     <<", muF12="<<rdai.m_muf12<<", muF22="<<rdai.m_muf22;
    return s;
  }

  struct ME_Weight_Info {
  public:
    mewgttype::code m_type;
    double m_B, m_VI, m_KP, m_RS;
    std::vector<double> m_wren,m_wfac;
    double m_x1, m_x2, m_y1, m_y2, m_mur2, m_muf2;
    int m_fl1, m_fl2;
    size_t m_oqcd, m_oew;
    Cluster_Sequence_Info m_clusseqinfo;
    std::vector<DADS_Info> m_dadsinfos;
    std::vector<RDA_Info> m_rdainfos;
    ME_Weight_Info():
      m_type(mewgttype::none),
      m_B(0.0), m_VI(0.0), m_KP(0.0), m_RS(0.0),
      m_wren(2,0.0), m_wfac(16,0.0),
      m_x1(1.0), m_x2(1.0), m_y1(1.0), m_y2(1.0),
      m_mur2(0.0), m_muf2(0.0), m_fl1(0), m_fl2(0),
      m_oqcd(0), m_oew(0), m_clusseqinfo() {}
    ME_Weight_Info& operator*=(const double &scal);
    void Reset();
  };

  std::ostream & operator<<(std::ostream & s,
                            const ATOOLS::ME_Weight_Info & mwi);
}

#endif

