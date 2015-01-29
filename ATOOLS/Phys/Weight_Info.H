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
    s<<" w = "<<wi.m_weight<<std::endl;
    return s;
  }

  struct mewgttype {
    enum code {
      none    = 0,
      muR     = 1,
      muF     = 2
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

  inline std::ostream & operator<<(std::ostream & s,
                                   const mewgttype::code & type)
  {
    if (type==mewgttype::none) s<<"none";
    if (type==mewgttype::muR)  s<<"muR";
    if (type==mewgttype::muF)  s<<"muF";
    return s;
  }

  struct DADS_Info {
  public:
    double   m_wgt, m_mur2;
    PDF_Info m_pdf;
    DADS_Info(const double &wgt,
              const int &fl1=0, const int &fl2=0,
              const double &x1=0.,const double &x2=0.,
              const double &xf1=0.,const double &xf2=0.,
              const double &mur2=0.,const double &muf2=0.):
    m_wgt(wgt), m_mur2(mur2),
    m_pdf(PDF_Info(fl1,fl2,x1,x2,muf2,muf2,xf1,xf2)) {}
  };

  inline std::ostream & operator<<(std::ostream & s,
                                   const ATOOLS::DADS_Info & dadsi)
  {
    s<<"DADS: wgt="<<dadsi.m_wgt<<", muR2="<<dadsi.m_mur2<<", "
     <<dadsi.m_pdf;
    return s;
  }

  struct ME_Weight_Info {
  public:
    mewgttype::code m_type;
    double m_B, m_VI, m_KP, m_RS;
    std::vector<double> m_wren,m_wfac;
    double m_x1, m_x2, m_y1, m_y2, m_mur2, m_muf2;
    std::vector<DADS_Info> m_dadsinfos;
    ME_Weight_Info():
      m_type(mewgttype::none),
      m_B(0.0), m_VI(0.0), m_KP(0.0), m_RS(0.0),
      m_wren(2,0.0), m_wfac(16,0.0),
      m_x1(1.0), m_x2(1.0), m_y1(1.0), m_y2(1.0),
      m_mur2(0.0), m_muf2(0.0) {}
    ME_Weight_Info& operator*=(const double &scal);
  };

  std::ostream & operator<<(std::ostream & s,
                            const ATOOLS::ME_Weight_Info & mwi);
}

#endif

