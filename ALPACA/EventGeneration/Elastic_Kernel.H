#ifndef ALPACA_Event_Generation_Elastic_Kernel_H
#define ALPACA_Event_Generation_Elastic_Kernel_H

#include "ALPACA/EventGeneration/Simple_MEs.H"
#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
#include "ALPACA/EventGeneration/Analysis_Handler.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Phys/Flavour.H"
#include <vector>

namespace ALPACA {

  class Elastic_Kernel {

  private:
    std::shared_ptr<Dynamic_Quantities_Handler> p_dyn_quant_handler;
    std::shared_ptr<Analysis_Handler>           p_analysis_handler;

    Simple_MEs                m_mes;
    std::vector<double>       mu2_dyn;

    bool                      m_only_gluons; //Option to initialize for only gluons (exclude gg->qqb option), set in Event_Generator
    std::pair<bool, double>   m_fixed_sigma; //Option to use a fixed cross-section, set in Event_Generator
    
    double SelectThat_qiqi_qiqi(const double shat);
    double SelectThat_qiqj_qiqj(const double shat);
    double SelectThat_qiqbi_qiqbi(const double shat);
    double SelectThat_qiqbi_qjqbj(const double shat);
    double SelectThat_qiqbi_gg(const double shat);
    double SelectThat_qig_qig(const double shat);
    double SelectThat_gg_gg(const double shat);
    double SelectThat_gg_qiqbi(const double shat);
        
  public:
    Elastic_Kernel(std::shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler, 
                   std::shared_ptr<Analysis_Handler> ptr_analysis_handler);
    ~Elastic_Kernel();

    double operator()(const ATOOLS::Flavour in1, const ATOOLS::Flavour in2, 
		                  const ATOOLS::Flavour out1, const ATOOLS::Flavour out2, const double shat);
    double FindScatteringXsec(const ATOOLS::Flavour in1, const ATOOLS::Flavour in2, const double shat);
    double SelectThat(const ATOOLS::Flavour in1, const ATOOLS::Flavour in2, const double shat, 
		                  ATOOLS::Flavour & out1, ATOOLS::Flavour & out2);
    void   SetMu2Dyn(const std::vector<double> mu2_temp);

    std::vector<double> GetMu2Dyn() const { return mu2_dyn; }
  };
}
#endif
