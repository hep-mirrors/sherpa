#ifndef SHRIMPS_Cross_Sections_Sigma_Inelastic_H
#define SHRIMPS_Cross_Sections_Sigma_Inelastic_H
#include "SHRiMPS/Cross_Sections/Sigma_Base.H"
#include <map>

namespace SHRIMPS {
  
  struct eikcomp{
    bool operator() (const Omega_ik * eik1, const Omega_ik * eik2) const
    {
      if(eik1->FF1()->Number() < eik2->FF1()->Number()) return true;
      else if(eik1->FF1()->Number() > eik2->FF1()->Number()) return false;
      else {
	if (eik1->FF2()->Number() < eik2->FF2()->Number()) return true;
	else return false;
      }
    }
  };

  class Sigma_Inelastic : public Sigma_Base {
    double m_Bmin, m_Bmax, m_deltaB;
    std::map<Omega_ik *,double,eikcomp>               m_xsecs;
    std::map<Omega_ik *,std::vector<double>,eikcomp>  m_grids;
  public:
    Sigma_Inelastic() {};
    Sigma_Inelastic(std::list<Omega_ik *> * eikonals) : 
      Sigma_Base(eikonals) { }
    ~Sigma_Inelastic();
    void   SetSigma(const double & sigma=-1);
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    void   FillGrid(const double & Bmin,const double & Bmax,
		    const double & deltaB,const double & sigma);
    double FixEikonalAndImpact(Omega_ik *&);
    void   TestInelasticCrossSection();
 };
}
#endif
