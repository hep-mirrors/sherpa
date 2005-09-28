#ifndef Three_Body_MEs_H
#define Three_Body_MEs_H

#include "HD_ME_Base.H"

namespace HADRONS {
  class P_GammaFF : public HD_ME_Base {
    int m_phot,m_f1,m_f2;
  public:
    P_GammaFF(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
  };

  class P_2PGamma : public HD_ME_Base {
    int m_phot,m_p1,m_p2;
  public:
    P_2PGamma(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
  };

  class P_P2Gamma : public HD_ME_Base {
    int    m_phot1, m_phot2, m_p;
    double m_mrho2, m_grho2;
  public:
    P_P2Gamma(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
  };

  class P_3P_DalitzDef : public HD_ME_Base {
  private:
    bool m_allpions, m_allsame;
    int  m_pi0, m_pip, m_pim;
  public:
    P_3P_DalitzDef(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
  };

//  class P_3P_Dalitz : public HD_ME_Base {
//  private:
//    double m_a1,m_a2,m_b1,m_b2,m_c;
//  public:
//    P_3P_Dalitz(int,ATOOLS::Flavour *);
//    double operator()(const ATOOLS::Vec4D *);
//    void   SetDalitzParameters(std::vector<double> &);
//  };
//
  class P_3P_Dalitz : public HD_ME_Base {

    int    m_pseudo_1, m_pseudo_2, m_pseudo_3;
    double m_mProd1, m_mProd2, m_mProd3, m_mDecayer;
    double m_const,m_liny,m_linx,m_linyphase,m_linxphase,
      m_quady,m_quadx,m_quadyphase,m_quadxphase,
      m_phaseliny,m_phaselinx,m_phasequady,m_phasequadx;
  public:
    P_3P_Dalitz(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
    void SetModelParameters( GeneralModel _md );
  };
}


#endif
