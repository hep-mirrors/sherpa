#ifndef Primordial_KPerp_H
#define Primordial_KPerp_H

#include "Particle_List.H"
#include "Blob_List.H"
#include "Poincare.H"

namespace SHERPA {

  class Primordial_KPerp {
  private:
    
    std::vector<ATOOLS::Vec3D>  *p_kperp[2];
    std::set<ATOOLS::Particle*> *p_filled;
    std::set<ATOOLS::Blob*>     *p_boosted;

    ATOOLS::Poincare m_oldcms, m_rotate, m_newcms;

    int m_current[2], m_scheme;
    
    double m_kperpmean[2], m_kperpsigma[2];
    double m_kperpmeans[2], m_kperpsigmas[2];

    double Lambda2(double sp,double sp1,double sp2); 

    bool FindConnected(ATOOLS::Particle *particle,ATOOLS::Particle *&connected,
		       bool forward,unsigned int catcher);
    bool BoostConnected(ATOOLS::Blob*,unsigned int catcher);

  public:

    // constructor
    Primordial_KPerp(std::string _m_path,std::string _m_file);
    
    // destructor
    ~Primordial_KPerp();

    // setters
    void SetKPerpMean(const double _m_kperpmean,unsigned int beam);
    void SetKPerpSigma(const double _m_kperpsigma,unsigned int beam);
    void SetScheme(const int _m_scheme);

    // getters
    double       KPerpMean(unsigned int beam) const;
    double       KPerpSigma(unsigned int beam) const;
    unsigned int Scheme() const;

    // member functions
    bool CreateKPerp(ATOOLS::Blob *blob1,ATOOLS::Blob *blob2);

    void FillKPerp(ATOOLS::Particle *particle,unsigned int beam);
    void FillKPerp(ATOOLS::Blob *blob);

  };// end of class Primordial_KPerp

}// end of namespace SHERPA

#include "Primordial_KPerp.IC"

#endif
