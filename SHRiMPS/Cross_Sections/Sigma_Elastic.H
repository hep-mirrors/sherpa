#ifndef SHRIMPS_Cross_Sections_Sigma_Elastic_H
#define SHRIMPS_Cross_Sections_Sigma_Elastic_H
#include "SHRiMPS/Cross_Sections/Sigma_Base.H"

namespace SHRIMPS {
  class Sigma_Elastic : public Sigma_Base {
    class dSigma_dt : public Function_Base {
    protected:
      Sigma_Elastic * p_sigma_el;
      double          m_Q;
    public:
      dSigma_dt(Sigma_Elastic * sigma_el) : p_sigma_el(sigma_el) {}
      void   SetQ(const double & Q) { m_Q = Q; }
      double operator()(double B);
    };

    double m_logdelta;
    std::map<double,double> m_intgrid, m_diffgrid;

    void   FillDiffQGrid(const double & Qmax);
    double FillIntQGridAndNormalize();
  public:
    Sigma_Elastic() : m_logdelta(20.) {};

    void   FillDifferentialGrids();
    double SelectPT2() const;

    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetDiffArgument(const double & B);

    void   Test();
  };  
}
#endif
