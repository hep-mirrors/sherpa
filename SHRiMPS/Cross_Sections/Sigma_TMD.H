#ifndef SHRIMPS_Cross_Sections_Sigma_TMD_H
#define SHRIMPS_Cross_Sections_Sigma_TMD_H
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "REMNANTS/Tools/Pseudo_TMD.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class Sigma_TMD {
    REMNANTS::Pseudo_TMD   * p_tmd[2];
    MODEL::Strong_Coupling * p_alphaS;
    double m_xsec, m_dxsmax;
    double m_Ymax, m_S, m_pt02, m_pt2max_norm,m_xmin;
    double m_pt2[2], m_tau, m_y, m_ymin, m_ymax, m_x[2], m_phi[2], m_xtmd[2];
    ATOOLS::Vec4D m_P[2],m_qvec[2], m_kvec;
    double m_accu, m_sigma, m_maxdsigma;
    bool                                       m_ana;
    std::map<std::string,ATOOLS::Histogram * > m_histos;
    const double MakePoint();
    const double MakeMEWeight();
    const double gg2g();
    const bool   Calculate();
    void  SelectFlavours(const bool & fixflavour);
  public:
    Sigma_TMD(const double & pt02);
    ~Sigma_TMD();
    
    void Initialise(REMNANTS::Pseudo_TMD * tmd[2]); 
    const double MakeEvent();
    const double xTMD(const size_t beam,
		      const double & x,const double & Q2,
		      const double & pt2,const double & Y,
		      const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_gluon));
    inline void SetAlphaS(MODEL::Strong_Coupling * as) { p_alphaS = as; }
    inline void SetTMD(const size_t & beam,
		       REMNANTS::Pseudo_TMD * tmd)     { p_tmd[beam] = tmd; }
    inline const ATOOLS::Vec4D & GetProp(const size_t & beam) const { return m_qvec[beam]; }
    inline const double        & GetX(const size_t & beam)    const { return m_x[beam]; } 
    inline const double        & GetQT2(const size_t & beam)  const { return m_pt2[beam]; }
    inline const double        & xTMD(const size_t & beam)    const { return m_xtmd[beam]; } 
    inline const ATOOLS::Vec4D & GetEmit()                    const { return m_kvec; }
  };
}
#endif
