#ifndef ALPACA_Event_Generation_Split_Handler_H
#define ALPACA_Event_Generation_Split_Handler_H

#include "ALPACA/EventGeneration/Colour_Handler.H"
#include "ALPACA/Tools/Parton.H"
#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ALPACA/EventGeneration/Inelastic_Kernel.H"
#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
#include "ALPACA/EventGeneration/Analysis_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ALPACA/Tools/Split_Variables.H"

#include <memory>
#include <vector>
#include <iostream>

namespace ALPACA {

  typedef std::pair<std::shared_ptr<Parton>,std::shared_ptr<Parton>>            partpair;
  typedef std::pair<std::pair<double, std::pair<int, double>>, ATOOLS::Vec4D>   processpair;
  typedef std::pair<std::vector<std::shared_ptr<Parton>>,processpair>           taubarpair;
  typedef std::pair<ATOOLS::Vec4D,ATOOLS::Vec4D>                                vec4dpair;
  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour>                            flavourpair;
  typedef std::pair<vec4dpair,flavourpair>                                      kinematicspair;
  typedef std::pair<std::vector<ATOOLS::Vec4D>,std::vector<double>>             energypair;
  typedef std::pair<energypair,flavourpair>                                     kinematicspair_2;
  typedef std::pair<partpair,kinematicspair_2>                                  splitpair;

  class Split_Handler {
  private :
    std::shared_ptr<std::vector<std::pair<double,taubarpair>>>              p_taubars;
    std::shared_ptr<std::vector<std::pair<double,splitpair>>>               p_tausplits;
    
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>   p_partons;
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>   p_removed_partons;

    std::shared_ptr<Dynamic_Quantities_Handler>           p_dyn_quant_handler;
    std::shared_ptr<Analysis_Handler>                     p_analysis_handler;
    std::shared_ptr<Inelastic_Kernel>                     p_inelastic_kernel;
    std::shared_ptr<Colour_Handler>                       p_colour_handler;

    std::shared_ptr<std::pair<double, int>>  p_tau_restart_list;

    
    double                     m_taumax;
    double                     m_p_min;
    double                     m_OE_mult_split;
    double                     m_taumax_scaling_limit;
    double                     m_tsample_min, m_tsample_max;

    int                        m_timekeeper;
    
    bool                       m_splitting_merging;
    bool                       m_only_gluons;
    bool                       m_include_bose_factors;
    bool                       m_elastic_scattering;

    std::pair<bool, double>    m_fixed_gamma;

    //To be removed using new classes instead
    std::shared_ptr<std::pair<std::shared_ptr<Parton>, kinematicspair>>  p_save_split_kinematics;

    //Misc. functions
    double          Abs3(ATOOLS::Vec4D v);
    
  public:
    Split_Handler(std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_partons,
                  std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_removed_partons,
                  std::shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler,
                  std::shared_ptr<Analysis_Handler> ptr_analysis_handler,
                  std::shared_ptr<Inelastic_Kernel> ptr_inelastic_kernel,
                  std::shared_ptr<std::pair<double, int>> ptr_tau_restart_list,
                  std::shared_ptr<std::pair<std::shared_ptr<Parton>, kinematicspair>>  ptr_save_split_kinematics,
                  std::shared_ptr<std::vector<std::pair<double,taubarpair>>> ptr_taubars,
                  std::shared_ptr<std::vector<std::pair<double,splitpair>>> ptr_tausplits,
                  std::shared_ptr<Colour_Handler> ptr_colour_handler);
    ~Split_Handler();

    //Main functions
    std::pair<double,splitpair>                      FindTauSplit(std::shared_ptr<Parton> part_in, double tau);
    std::vector<double>                              FindSplittingOE(std::shared_ptr<Parton> part_in, double T_star, double mg2); 
    std::pair<std::vector<ATOOLS::Vec4D>, double>    FindSplittingKinematics(std::shared_ptr<Parton> part_in, std::shared_ptr<Parton> part_recoil, double x, double mg2, double kT2);
    std::pair<double, partpair>                      DoSplittingKinematics(double tausplit, splitpair split_pair);
    std::pair<bool, splitpair>                       DoesSplit(double tausplit, splitpair split_pair);
  };
}

#endif
