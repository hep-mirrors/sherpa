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

  private:    
    double m_tmin, m_tmax, m_summed;
    size_t m_steps;
    double m_delta;
    std::vector<double> m_intgrid, m_diffgrid;

    void FillDiffQGrid();
    void FillIntQGridAndNormalize();
  public:
    Sigma_Elastic(); 
    void   FillGrids();
    double SelectT() const;

    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetDiffArgument(const double & B);
    const double & Summed() const { return m_summed; }
    const double & Tmin()   const { return m_tmin; }   
    const double & Tmax()   const { return m_tmax; }   
    const size_t & Steps()  const { return m_steps; }
    const std::vector<double> & GetDiffGrid() const { return m_diffgrid; }
    double Test();
    double GetXSvsT(double t);

  };  
}
#endif
