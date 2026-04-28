#ifndef BEAM_Ions_Ion_Base_H
#define BEAM_Ions_Ion_Base_H

#include "BEAM/Main/Beam_Base.H"
#include <list>


namespace BEAM {
  struct Ion_Constituent {
    ATOOLS::Flavour m_flav;
    ATOOLS::Vec4D   m_r, m_p;
    double          m_E, m_q2;
    int             m_Rform, m_Pform;
    
    Ion_Constituent(const ATOOLS::Flavour & flav,
		    const ATOOLS::Vec4D & r=ATOOLS::Vec4D(),
		    const ATOOLS::Vec4D & p=ATOOLS::Vec4D(),
		    const double E=0.) :
      m_flav(flav), m_r(r), m_p(p), m_E(E), m_q2(m_p.Abs2()) {} 
  };

  class Ion_Propagation;
  
  class Ion_Base : public Beam_Base {
  protected:
    size_t          m_A, m_Z, m_Np, m_Nn;
    double          m_R, m_rho0, m_kFermi, m_rmin2, m_dmin2;
    int             m_Rform, m_Pform;
    std::vector<Ion_Constituent *> m_constituents;

    Ion_Propagation * p_prop;
    
    inline const double kernel(const double & kappa) const {
      double f = m_beam.Mass();
      for (std::vector<Ion_Constituent *>::const_iterator icit=m_constituents.begin();
	   icit!=m_constituents.end();icit++) {
	f -= sqrt(ATOOLS::Vec3D((*icit)->m_p).Sqr()+
		  kappa*ATOOLS::sqr((*icit)->m_flav.Mass()));
      }
      return f;
    }
    inline const double dkernel(const double & kappa) const {
      double df = 0., m;
      for (std::vector<Ion_Constituent *>::const_iterator icit=m_constituents.begin();
	   icit!=m_constituents.end();icit++) {
	m   = (*icit)->m_flav.Mass();
	df -= m/sqrt(ATOOLS::Vec3D((*icit)->m_p).Sqr()+kappa*ATOOLS::sqr(m));
      }
      return df;
    }
    
    void Initialise();
    void MakeInitialPositions(bool withBlocking=true);
    ATOOLS::Vec4D MakeInitialPositionsHardSphere(bool withBlocking=true);
    ATOOLS::Vec4D MakeInitialPositionsGaussian(bool withBlocking=true);
    bool IsRBlocked(const ATOOLS::Flavour & flav,const ATOOLS::Vec4D & pos);
    bool IsPRBlocked(const size_t ref,const ATOOLS::Flavour & flav,
		     const ATOOLS::Vec4D & pos,const ATOOLS::Vec4D & mom);
    void MakeInitialMomenta(bool withBlocking=true);
    void FixEnergies();
    void PrintConstituents();
  public:
    Ion_Base(beamspectrum type,const ATOOLS::Flavour & beam,
	     const double energy,const double polarisation,
	     const int dir, int mode = 0);

    Beam_Base * Copy() {
      return new Ion_Base(m_type,m_beam,m_energy,m_polarisation,m_dir);
    };
    bool CalculateWeight(double, double) { return true; }
    void CheckRadius();
    
    inline const size_t & A()  const { return m_A; }
    inline const size_t & Z()  const { return m_Z; }
    inline const double & R()  const { return m_R; }
    inline const size_t & Np() const { return m_Np; }
    inline const size_t & Nn() const { return m_Nn; }
    inline std::vector<Ion_Constituent *> * GetConstituents() {
      return &m_constituents;
    }
  };
}

#endif
