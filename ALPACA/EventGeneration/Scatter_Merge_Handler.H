#ifndef ALPACA_Event_Generation_Scatter_Merge_Handler_H
#define ALPACA_Event_Generation_Scatter_Merge_Handler_H

#include "ALPACA/EventGeneration/Colour_Handler.H"
#include "ALPACA/Tools/Parton.H"
#include "ALPACA/EventGeneration/Elastic_Kernel.H"
#include "ALPACA/EventGeneration/Inelastic_Kernel.H"
#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
#include "ALPACA/EventGeneration/Analysis_Handler.H"
#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Org/CXXFLAGS.H"

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

  class Scatter_Merge_Handler {
  private :
    std::shared_ptr<std::vector<std::pair<double,taubarpair>>>              p_taubars;
    std::shared_ptr<std::vector<std::pair<double,splitpair>>>               p_tausplits;
    
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>   p_partons;
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>   p_removed_partons;

    std::shared_ptr<Colour_Handler>                       p_colour_handler;
    std::shared_ptr<Dynamic_Quantities_Handler>           p_dyn_quant_handler;
    std::shared_ptr<Analysis_Handler>                     p_analysis_handler;
    std::shared_ptr<Inelastic_Kernel>                     p_inelastic_kernel;
    std::shared_ptr<Elastic_Kernel>                       p_elastic_kernel;

    std::shared_ptr<std::pair<double, int>>  p_tau_restart_list;


    
    double                     m_taumax;
    double                     m_p_min;
    double                     m_OE_mult_scatter;
    double                     m_OE_mult_merge;
    double                     m_taumax_scaling_limit;

    bool                       m_only_gluons;
    bool                       m_include_bose_factors;
    bool                       m_elastic_scattering;
    bool                       m_splitting_merging;

    int                        m_timekeeper;
    
    std::pair<bool, double>    m_fixed_sigma;
    

    //To be removed using new classes instead
    std::shared_ptr<std::pair<std::shared_ptr<Parton>, kinematicspair>>   p_save_merging_kinematics;
    std::shared_ptr<std::pair<kinematicspair, std::pair<double,double>>>  p_save_scatter_kinematics;

    //Misc. functions
    double          Abs3(ATOOLS::Vec4D v);
    
  public:
    Scatter_Merge_Handler(std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_partons,
                          std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_removed_partons,
                          std::shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler,
                          std::shared_ptr<Analysis_Handler> ptr_analysis_handler,
                          std::shared_ptr<Inelastic_Kernel> ptr_inelastic_kernel,
                          std::shared_ptr<Elastic_Kernel> ptr_elastic_kernel,
                          std::shared_ptr<std::pair<double, int>> ptr_tau_restart_list,
                          std::shared_ptr<std::pair<std::shared_ptr<Parton>, kinematicspair>>  ptr_save_merging_kinematics,
                          std::shared_ptr<std::pair<kinematicspair, std::pair<double,double>>>  ptr_save_scatter_kinematics,
                          std::shared_ptr<std::vector<std::pair<double,taubarpair>>> ptr_taubars,
                          std::shared_ptr<std::vector<std::pair<double,splitpair>>> ptr_tausplits,
                          std::shared_ptr<Colour_Handler> ptr_colour_handler);

    ~Scatter_Merge_Handler();

    //Main functions
    std::pair<double,taubarpair>                 FindTauScatterMerge(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double tau);
    std::pair<int, std::pair<double, double>>    DoesScatterOrMerge(const double taubar, taubarpair taubar_pair);
    bool                                         FindScatteringKinematics(std::shared_ptr<Parton> part1, std::shared_ptr<Parton> part2, const double taubar, double mg2, double mf2);
    int                                          DoScatteringKinematics(std::shared_ptr<Parton> part1, std::shared_ptr<Parton> part2, const double taubar, ATOOLS::Vec4D shift);
    std::pair<bool, std::vector<double>>         FindMergingKinematics(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, std::shared_ptr<Parton> part_recoil);
    std::pair<bool,std::shared_ptr<Parton>>      DoMergingKinematics(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double taumerge, double t_form);

  };
}

#endif
