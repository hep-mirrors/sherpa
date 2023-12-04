#ifndef ALARIC__Shower__Cluster_H
#define ALARIC__Shower__Cluster_H

#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ALARIC++/Tools/Splitting.H"

namespace ALARIC {

  class Kernel;
  class Shower;

  struct CS_Parameters {
    double m_kt2, m_z, m_y, m_phi, m_x, m_wk, m_ws, m_q2, m_mu2;
    int m_mode, m_kin, m_col, m_kmode;
    Kernel *p_sf;
    size_t m_idi, m_idj, m_idk;
    ATOOLS::Vec4D m_pbt;
    CS_Parameters(const double &kt2,const double &z,
		  const double &y,const double &phi,
		  const double &x,const double &q2,
		  const int mode=-1,const int kin=0,const int kmode=0):
      m_kt2(kt2), m_z(z), m_y(y), m_phi(phi), m_x(x), m_wk(1.0), m_ws(kt2),
      m_q2(q2), m_mu2(kt2), m_mode(mode), m_kin(kin), m_col(0), m_kmode(kmode), p_sf(NULL),
      m_idi(0), m_idj(0), m_idk(0) {}
  };// end of struct CS_Parameters

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs);

  class Cluster_Definitions: public PDF::Cluster_Definitions_Base {
  private:

    const ATOOLS::Mass_Selector *p_ms;

    ATOOLS::Cluster_Leg *p_b;
    
    Shower *p_shower;
    int     m_mode, m_amode;

    double Lambda(const double &s,const double &sb,const double &sc) const;
    double Phi(ATOOLS::Vec4D pijt,ATOOLS::Vec4D pkt,
	       ATOOLS::Vec4D pi,const bool ii=false) const;
    ATOOLS::Flavour ProperFlav(const ATOOLS::Flavour &fl) const;

  public:
    
    Cluster_Definitions(Shower *const shower);

    PDF::Cluster_Param Cluster(const PDF::Cluster_Config &ca);

    Splitting KT2(const ATOOLS::Cluster_Amplitude &ampl,
		  int i,int j,int k,const ATOOLS::Flavour &mo,
		  const int id,const int type,const int mode,
		  double &ws,double &mu2);
    
    inline ATOOLS::Cluster_Leg *LegB() const { return p_b; }

    inline void SetAMode(const size_t mode) { m_amode=mode; }

  };// end of class Cluster_Definitions

}// end of namespace ALARIC

#endif
