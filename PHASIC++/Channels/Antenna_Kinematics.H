#ifndef PHASIC__Channels__Antenna_Kinematics_H
#define PHASIC__Channels__Antenna_Kinematics_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"

namespace PHASIC {

  struct Ant_Args {
    double m_y, m_z, m_phi;
    ATOOLS::Vec4D m_pi, m_pj, m_pk, m_pijt, m_Kt, m_nb, m_K, m_n;
    ATOOLS::Vec4D_Vector m_p;
    std::vector<int> m_b;
    int m_mode, m_stat, m_ij, m_k;
    inline Ant_Args(const double &y=0.0,const double &z=0.0,
		    const double &phi=0.0,const int mode=0,
		    const int stat=-1):
      m_y(y), m_z(z), m_phi(phi),
      m_mode(mode), m_stat(stat) {}
  };// end of struct Ant_Args

  int ClusterAntenna(Ant_Args &ffp,const size_t &i,const size_t &j,
		     const size_t &k,const double &mij2);

  int ConstructAntenna(Ant_Args &ffp,const size_t &ij,const size_t &k,
		       const double &mi2,const double &mj2);

}// end of namespace PHASIC

#endif
