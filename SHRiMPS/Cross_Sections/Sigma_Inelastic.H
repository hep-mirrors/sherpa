#ifndef SHRIMPS_Cross_Sections_Sigma_Inelastic_H
#define SHRIMPS_Cross_Sections_Sigma_Inelastic_H
#include "SHRiMPS/Cross_Sections/Sigma_Base.H"
#include <map>

namespace SHRIMPS {
  class Sigma_Inelastic : public Sigma_Base {
  private:
    std::map<Omega_ik *,double>                 m_xsecs; 
    std::map<Omega_ik *,std::vector<double> * > m_Bgrids; 

    const double FillBGrid(Omega_ik * eikonal);
  public:
    Sigma_Inelastic() {}
    ~Sigma_Inelastic();

    void       FillDifferentialGrids();
    Omega_ik * SelectEikonal();
    double     SelectB();
    
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double XSec(Omega_ik * eik=NULL);
    double Test();
 };
}
#endif
