#ifndef CSSHOWER_Main_CS_Cluster_Definitions_H
#define CSSHOWER_Main_CS_Cluster_Definitions_H

#include "ATOOLS/Phys/Cluster_Definitions_Base.H"

namespace CSSHOWER {

  class Splitting_Function_Base;
  class Shower;

  struct CS_Parameters {
    double m_kt2, m_z, m_y, m_phi, m_x, m_wk, m_ws;
    int m_mode;
    CS_Parameters(const double &kt2,const double &z,
		  const double &y,const double &phi,
		  const double &x,const int mode=-1):
      m_kt2(kt2), m_z(z), m_y(y), m_phi(phi), m_x(x),
      m_wk(1.0), m_ws(kt2), m_mode(mode) {}
  };// end of struct CS_Parameters

  class CS_Cluster_Definitions: public ATOOLS::Cluster_Definitions_Base {
  private:

    ATOOLS::Mass_Selector *p_ms;
    ATOOLS::Cluster_Leg   *p_b;
    
    Shower *p_shower;
    int     m_mode, m_kmode;

    double Lambda(const double &s,const double &sb,const double &sc) const;
    double Phi(ATOOLS::Vec4D pijt,ATOOLS::Vec4D pkt,
	       ATOOLS::Vec4D pi,const bool ii=false) const;
    void KernelWeight(const ATOOLS::Cluster_Leg *i,
		      const ATOOLS::Cluster_Leg *j,
		      const ATOOLS::Cluster_Leg *k,
		      const ATOOLS::Flavour &mo,CS_Parameters &cs) const;
    double GetX(const ATOOLS::Cluster_Leg *l,
		Splitting_Function_Base *const sf) const;
    ATOOLS::Flavour ProperFlav(const ATOOLS::Flavour &fl) const;

  public:
    
    CS_Cluster_Definitions(Shower *const shower,const int kmode=1);

    // CS cluster kperps
    ATOOLS::CParam KPerp2(const ATOOLS::Cluster_Amplitude &ampl,
		  int i,int j,int k,
		  const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		  ATOOLS::Mass_Selector *const ms=NULL);
    CS_Parameters KT2(const ATOOLS::Cluster_Amplitude *ampl,
		      const ATOOLS::Cluster_Leg *i,
		      const ATOOLS::Cluster_Leg *j,
		      const ATOOLS::Cluster_Leg *k,
		      const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
		      ATOOLS::Mass_Selector *const ms=NULL);

    CS_Parameters KT2_FF(const ATOOLS::Cluster_Leg *i,
			 const ATOOLS::Cluster_Leg *j,
			 const ATOOLS::Cluster_Leg *k,
			 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    CS_Parameters KT2_FI(const ATOOLS::Cluster_Leg *i,
			 const ATOOLS::Cluster_Leg *j,
			 const ATOOLS::Cluster_Leg *a,
			 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    CS_Parameters KT2_IF(const ATOOLS::Cluster_Leg *a,
			 const ATOOLS::Cluster_Leg *i,
			 const ATOOLS::Cluster_Leg *k,
			 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    CS_Parameters KT2_II(const ATOOLS::Cluster_Leg *a,
			 const ATOOLS::Cluster_Leg *i,
			 const ATOOLS::Cluster_Leg *b,
			 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));

    // CS cluster kinematics
    ATOOLS::Vec4D_Vector Combine(const ATOOLS::Cluster_Amplitude &ampl,
				 int i,int j,int k,
				 const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none),
				 ATOOLS::Mass_Selector *const ms=NULL);
    
    ATOOLS::Vec4D_Vector Combine_FF(const ATOOLS::Cluster_Amplitude &ampl,
				    int i,int j,int k,
				    const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    ATOOLS::Vec4D_Vector Combine_FI(const ATOOLS::Cluster_Amplitude &ampl,
				    int i,int j,int k,
				    const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    ATOOLS::Vec4D_Vector Combine_IF(const ATOOLS::Cluster_Amplitude &ampl,
				    int i,int j,int k,
				    const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    ATOOLS::Vec4D_Vector Combine_II(const ATOOLS::Cluster_Amplitude &ampl,
				    int i,int j,int k,
				    const ATOOLS::Flavour &mo=ATOOLS::Flavour(kf_none));
    
    //inline methods
    inline double Lambda(const double a, const double b, const double c) {
      return a*a+b*b+c*c-2.*a*b-2.*a*c-2.*b*c;
    }
  };// end of class CS_Cluster_Definitions

}// end of namespace CSSHOWER

#endif
