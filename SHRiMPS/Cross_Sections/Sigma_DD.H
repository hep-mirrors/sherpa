#ifndef SHRIMPS_Cross_Sections_Sigma_DD_H
#define SHRIMPS_Cross_Sections_Sigma_DD_H
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"
#include "SHRiMPS/Cross_Sections/Sigma_SD.H"

namespace SHRIMPS {
  class Sigma_DD : public Sigma_Base {
    class dSigma_dt_Kernel : public Function_Base {
    protected:
      Omega_ik * p_eikonal; 
      double     m_Q;
    public:
      dSigma_dt_Kernel() {}
      void   SetQ(const double & Q)         { m_Q = Q; }
      void   SetEikonal(Omega_ik * eikonal) { p_eikonal = eikonal; }
      double operator()(double B);
    };
    Sigma_Elastic     * p_sigma_el;
    Sigma_SD          * p_sigma_sd;


    double              m_Bmin, m_Bmax;
    double              m_Qmax;
    size_t              m_logQsteps;
    double              m_logdelta;
    std::vector<double> m_intgrid_DD, m_diffgrid_DD;
    const std::vector<double> *p_grid;
    double              m_sigma_DD;
    int                 m_test;

    void FillGrids();
  public:
    Sigma_DD() {};
    Sigma_DD(Sigma_Elastic * sigma_el,Sigma_SD * sigma_sd,const int & test=0);
    //double Calculate(const double & Bmin, const double & Bmax,
    //		     Omega_ik * eikonal=NULL) {
    // FillGrids();
    //  return m_sigma_SD1+m_sigma_SD2;
    // }
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetDiffArgument(const double & B);
    double PT2();
    void   PrintDifferentialElasticAndDiffXsec(const bool & onscreen=false,
					       std::string dirname=std::string("."));
  };  
}
#endif
