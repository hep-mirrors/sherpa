#ifndef ALPACA_Event_Generation_Inelastic_Kernel_H
#define ALPACA_Event_Generation_Inelastic_Kernel_H

#include "ALPACA/Tools/Parton.H"
#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
#include "ALPACA/EventGeneration/Analysis_Handler.H"
#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#include <memory>
#include <vector>
#include <iostream>

namespace ALPACA {

  class Inelastic_Kernel {
  private :
    std::shared_ptr<Dynamic_Quantities_Handler> p_dyn_quant_handler;
    std::shared_ptr<Analysis_Handler>           p_analysis_handler;

    std::pair<bool, double>   m_fixed_gamma;
    bool                      m_gaussian_kT2;
    double                    m_kT2_reg;
    double                    m_formation_time;
    
  public:
    Inelastic_Kernel(std::shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler, std::shared_ptr<Analysis_Handler> ptr_analysis_handler);
    ~Inelastic_Kernel();

    double                    FindGamma(double p, const ATOOLS::Flavour flavour_a, const ATOOLS::Flavour flavour_b, const ATOOLS::Flavour flavour_c, double x, double Tstar, double mg2, double mf2);
    double                    FindCs1s2s3(ATOOLS::Flavour s1, ATOOLS::Flavour s2, ATOOLS::Flavour s3);
    double                    QBH(double r);
    std::complex<double>      Li2(std::complex<double> z, double tol, int N_max);
    double                    FindFormationTime(double E_in, double T_star, double mg2, double kT2);
    double                    SamplekT2(double kT2_min, double kT2_max);
    double                    FindMergingXsec(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double T_star, double mg2, double mf2, double x, double kT2, double Ep, double Ek, double p_out_com, double kT2_max);
    
  };
}

#endif
