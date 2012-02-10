#ifndef SINIC_Cross_Sections_Sigma_Total_H
#define SINIC_Cross_Sections_Sigma_Total_H
#include "SINIC++/Cross_Sections/Sigma_Base.H"

namespace SINIC {
  class Sigma_Tot : public Sigma_Base {
  public:
    Sigma_Tot(std::list<Omega_ik *> * eikonals,
	      const double & originalY, const double & cutoffY,
	      const double & accu) : 
      Sigma_Base(eikonals) { }
    double GetValue(const double & B);
    double GetCombinedValue(const double & B);
    void TestTotalCrossSection();
  };

  class Elastic_Slope : public Sigma_Base {
    double m_stot;
  public:
    Elastic_Slope(std::list<Omega_ik *> * eikonals,
	      const double & originalY, const double & cutoffY,
		  const double & accu,const double & stot) : 
      Sigma_Base(eikonals), m_stot(stot) { }
    double GetValue(const double & B);
    double GetCombinedValue(const double & B);
  };
}
#endif
