#ifndef SHERPA_SoftPhysics_Primordial_KPerp_H
#define SHERPA_SoftPhysics_Primordial_KPerp_H

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Poincare.H"

namespace PDF    { class Remnant_Base; }

namespace SHERPA {

  class Primordial_KPerp {
  private:
    
    std::vector<ATOOLS::Vec3D>  *p_kperp[2];
    std::set<ATOOLS::Particle*> *p_filled;
    std::set<ATOOLS::Blob*>     *p_boosted;

    ATOOLS::Poincare m_oldcms, m_rotate, m_newcms;

    int m_current[2], m_scheme, m_fill;
    
    double m_kperpmean[2], m_kperpsigma[2];
    double m_kperpmeans[2], m_kperpsigmas[2];

    PDF::Remnant_Base *p_remnants[2];

    double Lambda2(double sp,double sp1,double sp2); 

    bool FindConnected(ATOOLS::Particle *particle,ATOOLS::Particle *&connected,
		       bool forward,unsigned int catcher);
    bool BoostConnected(ATOOLS::Blob*,unsigned int catcher);

    void FillKPerp(ATOOLS::Particle *particle,unsigned int beam);

  public:

    // constructor
    Primordial_KPerp(std::string _m_path,std::string _m_file);
    
    // destructor
    ~Primordial_KPerp();

    // member functions
    bool CreateKPerp(ATOOLS::Blob *blob1,ATOOLS::Blob *blob2);

    void FillKPerp(ATOOLS::Blob *blob);

    // inline functions
    inline void SetKPerpMean(const double _m_kperpmean,unsigned int beam)
    { if (beam<2) m_kperpmean[beam]=ATOOLS::dabs(_m_kperpmean); }
    inline void SetKPerpSigma(const double _m_kperpsigma,unsigned int beam)
    { if (beam<2) m_kperpsigma[beam]=ATOOLS::dabs(_m_kperpsigma); }
    inline void SetRemnant(PDF::Remnant_Base *const rem,const size_t beam)
    { if (beam<2) p_remnants[beam]=rem; }
    inline void SetScheme(const int _m_scheme)
    { m_scheme=_m_scheme; }

    inline double KPerpMean(unsigned int beam) const
    { if (beam<2) return m_kperpmean[beam]; else return -1.0; }
    inline double KPerpSigma(unsigned int beam) const
    { if (beam<2) return m_kperpsigma[beam]; else return -1.0; }
    inline unsigned int Scheme() const
    { return m_scheme; }

  };// end of class Primordial_KPerp

}// end of namespace SHERPA

#endif
