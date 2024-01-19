#ifndef SHRIMPS_Cross_Sections_Sigma_D_H
#define SHRIMPS_Cross_Sections_Sigma_D_H
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"

namespace SHRIMPS {
  class Sigma_D : public Sigma_Base {
    class D_Term : public Function_Base {
    protected:
      Omega_ik * p_eikonal;
      double     m_Q;
    public:
      D_Term() {}
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

    //double GetCombinedValueEL(const double & B);
  public:
    Sigma_D();

    double SelectT(const size_t & mode) const;
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetValuePerChannel(const int i, const int j, const double & B);
    double GetCombinedValueSD0(const double & B);
    double GetCombinedValueSD1(const double & B);
    double GetCombinedValueDD(const double & B);
    double GetCombinedValueEL(const double & B);
    double GetXSec(const size_t diff) { return m_summed[diff]; }
    double GetXSvsT(const size_t diff, double t);
    void   FillGrids(Sigma_Elastic * sigma_el);
  };  
}
#endif
