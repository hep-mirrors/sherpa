#ifndef PHASIC__Channels__CSS_Kinematics_H
#define PHASIC__Channels__CSS_Kinematics_H

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  struct Kin_Args {
    double m_y, m_z, m_phi, m_mk2, m_x;
    ATOOLS::Vec4D m_pi, m_pj, m_pk, m_nb;
    ATOOLS::Poincare_Sequence m_lam;
    int m_mode, m_stat;
    static double s_uxeps;
    inline Kin_Args(const double &y=0.0,const double &z=0.0,
		    const double &phi=0.0,const int mode=0,
		    const int stat=-1):
      m_y(y), m_z(z), m_phi(phi), m_mk2(-1.0), m_x(0.0),
      m_mode(mode), m_stat(stat) {}
  };// end of struct Kin_Args

  Kin_Args ClusterFFDipole(const double &mi2,const double &mj2,
			   const double &mij2,const double &mk2,
			   const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
			   const ATOOLS::Vec4D &pk,const int mode=0);
  int ConstructFFDipole(const double &mi2,const double &mj2,
			const double &mij2,const double &mk2,
			const ATOOLS::Vec4D &pij,const ATOOLS::Vec4D &pk,
			Kin_Args &ffp);

  Kin_Args ClusterFIDipole(const double &mi2,const double &mj2,
			   const double &mij2,const double &ma2,
			   const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
			   const ATOOLS::Vec4D &pa,const int mode=0);
  int ConstructFIDipole(const double &mi2,const double &mj2,
			const double &mij2,const double &ma2,
			const ATOOLS::Vec4D &pij,const ATOOLS::Vec4D &pa,
			Kin_Args &fip);

  Kin_Args ClusterIFDipole(const double &ma2,const double &mj2,
			   const double &maj2,const double &mk2,
			   const double &mb2,const ATOOLS::Vec4D &pa,
			   const ATOOLS::Vec4D &pj,const ATOOLS::Vec4D &pk,
			   const ATOOLS::Vec4D &pb,const int mode=0);
  int ConstructIFDipole(const double &ma2,const double &mj2,
			const double &maj2,const double &mk2,
			const double &mb2,const ATOOLS::Vec4D &paj,
			const ATOOLS::Vec4D &pk,const ATOOLS::Vec4D &pb,
			Kin_Args &ifp);

  Kin_Args ClusterIIDipole(const double &ma2,const double &mi2,
			   const double &mai2,const double &mb2,
			   const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pi,
			   const ATOOLS::Vec4D &pb,const int mode=0);
  int ConstructIIDipole(const double &ma2,const double &mi2,
			const double &mai2,const double &mb2,
			const ATOOLS::Vec4D &pa,const ATOOLS::Vec4D &pb,
			Kin_Args &iip);

}// end of namespace PHASIC

#endif
