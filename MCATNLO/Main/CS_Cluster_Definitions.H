#ifndef MCATNLO_Main_CS_Cluster_Definitions_H
#define MCATNLO_Main_CS_Cluster_Definitions_H

#include "PDF/Main/Cluster_Definitions_Base.H"

namespace MCATNLO {

  class Splitting_Function_Base;
  class Shower;

  struct CS_Parameters {
    double m_kt2, m_z, m_y, m_phi, m_x, m_wk, m_ws, m_q2, m_mu2;
    int m_mode, m_kin;
    Splitting_Function_Base *p_sf;
    size_t m_idi, m_idj, m_idk;
    ATOOLS::Vec4D m_pbt;
    CS_Parameters(const double &kt2,const double &z,
		  const double &y,const double &phi,
		  const double &x,const double &q2,
		  const int mode=-1,const int kin=0):
      m_kt2(kt2), m_z(z), m_y(y), m_phi(phi), m_x(x), m_wk(1.0), m_ws(kt2),
      m_q2(q2), m_mu2(kt2), m_mode(mode), m_kin(kin), p_sf(NULL),
      m_idi(0), m_idj(0), m_idk(0) {}
  };// end of struct CS_Parameters

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs);

  class CS_Cluster_Definitions: public PDF::Cluster_Definitions_Base {
  private:

    ATOOLS::Mass_Selector *p_ms;
    ATOOLS::Cluster_Leg   *p_b;
    
    Shower *p_shower;
    int     m_mode, m_kmode, m_amode;

    double Lambda(const double &s,const double &sb,const double &sc) const;
    double Phi(ATOOLS::Vec4D pijt,ATOOLS::Vec4D pkt,
	       ATOOLS::Vec4D pi,const bool ii=false) const;
    void KernelWeight(const ATOOLS::Cluster_Leg *i,
		      const ATOOLS::Cluster_Leg *j,
		      const ATOOLS::Cluster_Leg *k,
		      const ATOOLS::Flavour &mo,CS_Parameters &cs) const;
    ATOOLS::Flavour ProperFlav(const ATOOLS::Flavour &fl) const;

  public:
    
    CS_Cluster_Definitions(Shower *const shower,const int kmode=1);

    // CS cluster kperps
    double GetX(const ATOOLS::Cluster_Leg *l,
		Splitting_Function_Base *const sf) const;
    PDF::CParam KPerp2(const ATOOLS::Cluster_Amplitude &ampl,
		  int i,int j,int k,
		  const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		       ATOOLS::Mass_Selector *const ms=NULL,const int kin=-1);
    CS_Parameters KT2(const ATOOLS::Cluster_Amplitude *ampl,
		      const ATOOLS::Cluster_Leg *i,
		      const ATOOLS::Cluster_Leg *j,
		      const ATOOLS::Cluster_Leg *k,
		      const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		      ATOOLS::Mass_Selector *const ms=NULL,const int ikin=1);

    // CS cluster kinematics
    ATOOLS::Vec4D_Vector Combine(const ATOOLS::Cluster_Amplitude &ampl,
				 int i,int j,int k,
				 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
				 ATOOLS::Mass_Selector *const ms=NULL,
				 const int kin=0);
    
    //inline methods
    inline ATOOLS::Cluster_Leg *LegB() const { return p_b; }

    inline void SetAMode(const size_t mode) { m_amode=mode; }

  };// end of class CS_Cluster_Definitions

}// end of namespace MCATNLO

#endif
