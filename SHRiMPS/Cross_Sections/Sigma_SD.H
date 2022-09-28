#ifndef SHRIMPS_Cross_Sections_Sigma_SD_H
#define SHRIMPS_Cross_Sections_Sigma_SD_H
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"

namespace SHRIMPS {
  class Sigma_SD : public Sigma_Base {
    class SD_Term : public Function_Base {
    protected:
      Omega_ik * p_eikonal;
      double     m_Q;
    public:
      SD_Term() {}
      void   SetEikonal(Omega_ik * eikonal) { p_eikonal = eikonal; }
      void   SetQ(const double & Q)         { m_Q = Q; }
      double operator()(double B);
    };

    double m_tmin, m_tmax, m_summed[3];
    size_t m_steps;
    double m_delta;
    std::vector<std::vector<std::vector<double> > > m_tgrids;
    std::vector<double> m_diffgrids[3], m_intgrids[3];
    
    void FillTGrids();
    void CombineTGrids(const size_t diff);
    void CreateIntGrids(const size_t diff,Sigma_Elastic * sigma_el);
  public:
    Sigma_SD();

    double SelectT(const size_t & mode) const;
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetXSec(const size_t diff) { return m_summed[diff]; }
    void   FillGrids(Sigma_Elastic * sigma_el);
  };  
}
#endif
