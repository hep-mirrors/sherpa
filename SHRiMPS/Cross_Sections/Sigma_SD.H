#ifndef SHRIMPS_Cross_Sections_Sigma_SD_H
#define SHRIMPS_Cross_Sections_Sigma_SD_H
#include "SHRiMPS/Cross_Sections/Sigma_Elastic.H"

namespace SHRIMPS {
  class Sigma_SD : public Sigma_Base {
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


    double              m_Bmin, m_Bmax;
    double              m_Qmax;
    size_t              m_logQsteps;
    double              m_logdelta;
    std::vector<double> m_intgrid_SD1, m_intgrid_SD2, m_diffgrid_SD1, m_diffgrid_SD2;
    const std::vector<double> *p_grid;
    double              m_sigma_SD1,   m_sigma_SD2;
    int                 m_test;

    void FillGrids();
  public:
    Sigma_SD() {};
    Sigma_SD(Sigma_Elastic * sigma_el,const int & test=0);
    //double Calculate(const double & Bmin, const double & Bmax,
    //		     Omega_ik * eikonal=NULL) {
    // FillGrids();
    //  return m_sigma_SD1+m_sigma_SD2;
    // }
    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double GetDiffArgument(const double & B);
    double PT(bool & mode);
    void   PrintDifferentialElasticAndSDXsec(const bool & onscreen=false,
					     std::string dirname=std::string("."));
    void   TestElasticCrossSection();
    const std::vector<double> * Grid1() const { return &m_intgrid_SD1; }
    const std::vector<double> * Grid2() const { return &m_intgrid_SD2; }
  };  
}
#endif
