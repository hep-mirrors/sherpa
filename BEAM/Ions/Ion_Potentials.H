#ifndef BEAM_Main_Ion_Potentials_H
#define BEAM_Main_Ion_Potentials_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <vector>

namespace BEAM {
  struct Ion_Constituent;
  class  Ion_Base;
  
  class Ion_Potentials {
  private:
    double m_A;
    double m_rho0, m_sigma, m_sigma2, m_norm, m_rhonorm;
    double m_V2loc, m_V3loc, m_maxRloc, m_nu, m_U2pref, m_U3pref;
    double m_aYuk, m_VYuk, m_maxRYuk, m_VYukpref, m_UYukpref;
    std::vector<std::vector<std::vector<ATOOLS::Vec3D> > > * p_deltarvec;
    std::vector<std::vector<std::vector<double> > >        * p_deltar, * p_deltar2norm;
    std::vector<double>                                      m_tilderho, m_rho, m_U;
  public:
    Ion_Potentials(size_t A);

    void UpdateDensities();
    void UpdatePotentials();

    const double LocalPotential2(const size_t & i);
    const double LocalPotential3(const size_t & i);
    const double YukawaPotential(const size_t & i);
    ATOOLS::Vec3D NablaLocalPotential2(const size_t & level,const size_t & i); 
    ATOOLS::Vec3D NablaLocalPotential3(const size_t & level,const size_t & i); 
    ATOOLS::Vec3D NablaYukawa(const size_t & level,const size_t & i); 
    ATOOLS::Vec3D Nabla(const size_t & level,const size_t & i); 

    
    inline const double & U(size_t i)        const { return m_U[i]; }
    inline const double & Rho(size_t i)      const { return m_rho[i]; }
    inline const double & RhoTilde(size_t i) const { return m_tilderho[i]; }
    inline void
    SetDistances(std::vector<std::vector<std::vector<ATOOLS::Vec3D> > > * deltarvec,
		 std::vector<std::vector<std::vector<double> > > * deltar,
		 std::vector<std::vector<std::vector<double> > > * deltar2norm) {
      p_deltarvec   = deltarvec;
      p_deltar      = deltar;
      p_deltar2norm = deltar2norm;
    }
  };
}

#endif
  
