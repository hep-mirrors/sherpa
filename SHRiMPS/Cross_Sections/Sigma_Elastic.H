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

    double              m_Bmin, m_Bmax;
    double              m_Qmax;
    size_t              m_logQsteps;
    double              m_logdelta;
    std::vector<double> m_intgrid, m_diffgrid;
    int                 m_test;

    void FillGrid();
  public:
    Sigma_Elastic() {};
    Sigma_Elastic(std::list<Omega_ik *> * eikonals,const double & energy,
		  const int & test=0);
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetDiffArgument(const double & B);
    double PT2() const;
    void   PrintDifferentialelasticXsec(const bool & onscreen=false, const bool & tuning=false, 
					std::string dirname=std::string("."));
    void   TestElasticCrossSection();

    const double & Bmin() const  { return m_Bmin; }
    const double & Bmax() const  { return m_Bmax; }
    const double & Qmax() const  { return m_Qmax; }
    const size_t & Steps() const { return m_logQsteps; }
    const double & Delta() const { return m_logdelta; }
    const std::vector<double> * Grid() const { return &m_intgrid; }
  };  
}
#endif
