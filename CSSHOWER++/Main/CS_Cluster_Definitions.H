#ifndef CSSHOWER_Main_CS_Cluster_Definitions_H
#define CSSHOWER_Main_CS_Cluster_Definitions_H

#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Process/Subprocess_Info.H"

namespace CSSHOWER {

  class Splitting_Function_Base;
  class Shower;

  struct CS_Parameters {
    double m_kt2, m_z, m_y, m_phi, m_x, m_wk, m_ws, m_q2, m_mu2;
    int m_mode, m_kin, m_col, m_kmode;
    Splitting_Function_Base *p_sf;
    ATOOLS::Vec4D m_pk;
    CS_Parameters(const double &kt2,const double &z,
		  const double &y,const double &phi,
		  const double &x,const double &q2,
		  const int mode,const int kin,const int kmode):
      m_kt2(kt2), m_z(z), m_y(y), m_phi(phi), m_x(x), m_wk(1.0), m_ws(kt2),
      m_q2(q2), m_mu2(kt2), m_mode(mode), m_kin(kin), m_col(0), m_kmode(kmode), p_sf(NULL) {}
  };// end of struct CS_Parameters

  std::ostream &operator<<(std::ostream &str,const CS_Parameters &cs);

  class CS_Cluster_Definitions: public PDF::Cluster_Definitions_Base {
  private:

    ATOOLS::Mass_Selector *p_ms;
    ATOOLS::Cluster_Leg   *p_b;
    
    Shower *p_shower;
    int     m_kmode, m_meweight, m_pdfcheck;

    Splitting_Function_Base *GetSF(const ATOOLS::Cluster_Leg *i,
				   const ATOOLS::Cluster_Leg *j,
				   const ATOOLS::Cluster_Leg *k,
				   const ATOOLS::Flavour &mo,CS_Parameters &cs) const;
    void KernelWeight(const ATOOLS::Cluster_Leg *i,
		      const ATOOLS::Cluster_Leg *j,
		      const ATOOLS::Cluster_Leg *k,
		      const ATOOLS::Flavour &mo,CS_Parameters &cs,
		      const int kmode) const;
    ATOOLS::Flavour ProperFlav(const ATOOLS::Flavour &fl) const;

    bool CheckColors(const ATOOLS::Cluster_Leg *li,
		     const ATOOLS::Cluster_Leg *lj,
		     const ATOOLS::Cluster_Leg *lk,
		     const ATOOLS::Flavour &mo) const;

    double Differential(ATOOLS::Cluster_Amplitude *const ampl,
			const int mode) const;

  public:
    
    CS_Cluster_Definitions(Shower *const shower,const int kmode,
			   const int meweight,const int pdfcheck);

    // CS cluster kperps
    double GetX(const ATOOLS::Cluster_Leg *l,
		Splitting_Function_Base *const sf) const;
    PDF::CParam KPerp2(const ATOOLS::Cluster_Amplitude &ampl,
		  int i,int j,int k,
		  const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		       ATOOLS::Mass_Selector *const ms=NULL,
		       const int kin=-1,const int mode=0);
    CS_Parameters KT2(const ATOOLS::Cluster_Amplitude *ampl,
		      const ATOOLS::Cluster_Leg *i,
		      const ATOOLS::Cluster_Leg *j,
		      const ATOOLS::Cluster_Leg *k,
		      const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		      ATOOLS::Mass_Selector *const ms=NULL,
		      const int kin=-1,const int mode=0,const int force=0);

    // CS cluster kinematics
    ATOOLS::Vec4D_Vector Combine(const ATOOLS::Cluster_Amplitude &ampl,
				 int i,int j,int k,
				 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
				 ATOOLS::Mass_Selector *const ms=NULL,
				 const int kin=0,const int mode=0);
    
  };// end of class CS_Cluster_Definitions

}// end of namespace CSSHOWER

#endif
