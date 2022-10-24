#ifndef PHASIC__Channels__Antenna_Kinematics_H
#define PHASIC__Channels__Antenna_Kinematics_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  struct Ant_Args {
    double m_y, m_z, m_phi;
    ATOOLS::Vec4D m_pi, m_pj, m_pk;
    ATOOLS::Vec4D_Vector m_p;
    std::vector<int> m_b;
    int m_mode, m_stat, m_ij, m_k;
    inline Ant_Args(const double &y=0.0,const double &z=0.0,
		    const double &phi=0.0,const int mode=0,
		    const int stat=-1):
      m_y(y), m_z(z), m_phi(phi),
      m_mode(mode), m_stat(stat) {}
  };// end of struct Ant_Args

  int ConstructFSAntenna(const double &mi2,const double &mj2,const double &mk2,
			 const ATOOLS::Vec4D &pij,const ATOOLS::Vec4D &pk,
			 Ant_Args &ffp);
  int ConstructISAntenna(const double &mi2,const double &mj2,const double &mk2,
			 const ATOOLS::Vec4D &pij,const ATOOLS::Vec4D &pk,
			 Ant_Args &ffp);

}// end of namespace PHASIC

#endif
