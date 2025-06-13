#ifndef ALPACA_EventGeneration_Event_Generator_H
#define ALPACA_EventGeneration_Event_Generator_H

#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ALPACA/EventGeneration/Colour_Handler.H"
#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
#include "ALPACA/EventGeneration/Inelastic_Kernel.H"
#include "ALPACA/EventGeneration/Elastic_Kernel.H"
#include "ALPACA/EventGeneration/Scatter_Merge_Handler.H"
#include "ALPACA/EventGeneration/Split_Handler.H"
#include "ALPACA/EventGeneration/Analysis_Handler.H"
#include "ALPACA/Tools/Parton.H"
#include "ALPACA/Tools/Split_Variables.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include <memory>
#include <vector>
#include <map>
#include <cmath>
#include <complex>
#include <stdio.h>


namespace ALPACA {

  typedef std::pair<double,double>                                                       distpair; 
  typedef std::pair<std::shared_ptr<Parton>,std::shared_ptr<Parton>>                     partpair;
  typedef std::pair<std::pair<double, std::pair<int, double>>, ATOOLS::Vec4D>            processpair;
  typedef std::pair<std::vector<std::shared_ptr<Parton>>,processpair>                    taubarpair;
  typedef std::pair<ATOOLS::Vec4D,ATOOLS::Vec4D>                                         vec4dpair;
  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour>                                     flavourpair;
  typedef std::pair<vec4dpair,flavourpair>                                               kinematicspair;
  typedef std::pair<std::vector<ATOOLS::Vec4D>,std::vector<double>>                      energypair;
  typedef std::pair<energypair,flavourpair>                                              kinematicspair_2;
  typedef std::pair<partpair,kinematicspair_2>                                           splitpair;

  class Event_Generator {
  private:
    std::shared_ptr<std::vector<std::pair<double,taubarpair>>>              p_taubars; //List of tau's where each particle pair has their closest approach, used to check 2->2 scatters
    std::shared_ptr<std::vector<std::pair<double,splitpair>>>               p_tausplit; //List of tau's where each particle has a possibility to split
    
    std::shared_ptr<std::vector<std::pair<double,std::shared_ptr<Split_Variables>>>>  p_tausplit_v2; //List of tau's where each particle has a possibility to split
    

    std::shared_ptr<Dynamic_Quantities_Handler>            p_dyn_quant_handler;
    std::shared_ptr<Inelastic_Kernel>                      p_inelastic_kernel;
    std::shared_ptr<Elastic_Kernel>                        p_elastic_kernel;
    std::shared_ptr<Scatter_Merge_Handler>                 p_scatter_merge_handler;
    std::shared_ptr<Split_Handler>                         p_split_handler;
    std::shared_ptr<Analysis_Handler>                      p_analysis_handler;

    double                                                 m_tau; //Current tau of the event
    double                                                 m_taumax; //Initial parameters of the run, are controlled by the runcard (through HIParams)
    double                                                 m_p_min; //Minimum value for p used in 1->2 splitting and 2->1 merging, 2->2 scattering and initial thermal sampling
    double                                                 m_taumax_scaling_limit; //TEMPORARY: Set upperlimit for how long simulation can go in tau, even though particles exists with t<taumax
    double                                                 m_extra_tau; //TEMPORARY: keep track of how many extra iterations of m_tau has been checked   
    double                                                 m_tsample_min, m_tsample_max; //Time interval for which t_avg to sample dN/dt statistics for 2->2, 1->2 and 2->1
    
    std::shared_ptr<Colour_Handler>                        p_colour_handler;
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>    p_partons; //List of all the current particles in the run
    std::shared_ptr<std::list<std::shared_ptr<Parton>>>    p_removed_partons; //List of all the removed particles (from splitting or merging) in the run
    std::shared_ptr<std::list<std::pair<int,int>>>         p_flow_backtrack; //List of colours to keep track of where the colour partner is in a blob not converted and sent to ALPACA
    
    bool                                                   m_only_gluons; //Toggles runs to be without quarks, only gluons
    bool                                                   m_no_process_end; //If true, no process left before taumax, will end the run
    bool                                                   m_do_kinematics; //Toggles if dynamics are updated or scatter/split/merge is just for bookkeeping
    bool                                                   m_show_progress_bar; //Toggle if progress bar is shown during run or not
    
    std::shared_ptr<bool>                                  p_add_histo;    
    std::shared_ptr<std::pair<double, int>>                p_tau_restart_list;
    
    std::shared_ptr<std::pair<kinematicspair, std::pair<double,double>>>  p_save_scatter_kinematics; //[[P_out_1, P_out_2], [Flav_out_1, Flav_out_2]]
    std::shared_ptr<std::pair<std::shared_ptr<Parton>, kinematicspair>>   p_save_split_kinematics; //[part_recoil, [[P_out,P_recoil] , [Flav_out, Flav_out]]]
    std::shared_ptr<std::pair<std::shared_ptr<Parton>, kinematicspair>>   p_save_merging_kinematics;
    
    double          m_test_double; //Variable only used for testing purposes
    bool            m_test_bool;   //Variable only used for testing purposes



    //Functions related to the main loop of the event generator
    void                    Reset();
    void                    FillTauTables();
    void                    ResetScatterMergeSplitList();
    double                  EvolveInTau();
    void                    UpdateTauTables(std::shared_ptr<Parton> part1, int b_type);

    //Other functions
    void                    SetTimekeeperLabFrame();
    void                    ProgressBar(float progress);

  public:
    Event_Generator(std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_parton_list, 
                    std::shared_ptr<std::list<std::pair<int,int>>> ptr_flow_backtrack);
    ~Event_Generator();

    int  GenerateEvent();

    static bool sortTaubars(std::pair<double,taubarpair> i, std::pair<double,taubarpair> j) { return (i.first<j.first); } //Function needed to order m_taubars with particle pairs
    static bool sortTausplits(std::pair<double,splitpair> i, std::pair<double,splitpair> j) { return (i.first<j.first); } //Function needed to order m_taubounds with particle pairs

    void CleanUp(const size_t & mode=0);
  };
}
#endif

