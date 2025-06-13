#ifndef ALPACA_Event_Generation_Analysis_Handler_H
#define ALPACA_Event_Generation_Analysis_Handler_H

#include "ALPACA/Tools/Parton.H"
#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ALPACA/EventGeneration/Dynamic_Quantities_Handler.H"
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

  class Analysis_Handler {
  private :
    //General variables
    std::shared_ptr<std::list<std::shared_ptr<Parton>>> p_partons;
    std::shared_ptr<std::list<std::shared_ptr<Parton>>> p_removed_partons;

    std::shared_ptr<std::map<std::string, std::shared_ptr<ATOOLS::Histogram>>>    p_histomap;
    std::shared_ptr<std::map<std::string, std::shared_ptr<ATOOLS::Histogram_2D>>> p_histomap2D;

    std::shared_ptr<std::vector<std::pair<double,taubarpair>>>              p_taubars;
    std::shared_ptr<std::vector<std::pair<double,splitpair>>>               p_tausplits;

    std::shared_ptr<bool>                    p_add_histo;
    std::shared_ptr<std::pair<double, int>>  p_tau_restart_list;

    std::string             m_seed_str;

    double                  m_tsample_min, m_tsample_max;
    double                  m_p_min;

    bool                    m_show_event_information;
    

    //Event by event variables
    int                     m_N_scatter_event;
    int                     m_N_split_event;
    int                     m_N_merge_event;
    int                     m_N_scatter_process_event[7]; //0:gg->gg, 
                                                          //1:qqb->gg       2:gq->gq,
                                                          //3:qiqi->qiqi,   4:qiqib->qiqib,
                                                          //5:qiqib->qjqjb, 6:qiqj->qiqj
    int                     m_N_split_process_event[3]; //0:g->gg, 1:g->qqb, 2:q->gq
    int                     m_N_merge_process_event[3]; //0:gg->g, 1:qqb->g, 2:gq->q
    
    std::vector<double>     m_N_scatter_merge_overdraw_event;
    std::vector<double>     m_N_split_overdraw_event;
    
    
    //Across all event variables
    double                                                 m_v2_w_total;
    double                                                 m_v2_total;
    std::vector<double>                                    m_PL_total;
    std::vector<double>                                    m_PT_total;
    std::vector<double>                                    m_PLPT_total;



  public:
    Analysis_Handler(std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_partons, 
                     std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_removed_partons,
                     std::shared_ptr<std::pair<double, int>> ptr_tau_restart_list,
                     std::shared_ptr<bool> ptr_add_histo,
                     std::shared_ptr<std::vector<std::pair<double,taubarpair>>> ptr_taubars,
                     std::shared_ptr<std::vector<std::pair<double,splitpair>>> ptr_tausplits);
    ~Analysis_Handler();

    //Main functions
    void      PrintEventInfo(bool pre_event, double tau);
    void      PrePostSampling(bool pre_event, double tau, std::shared_ptr<Dynamic_Quantities_Handler> dyn_quant_handler);

    void      Histo(std::string histo_name, double histo_val);
    void      Histo2D(std::string histo_name, double histo_val_1, double histo_val_2);

    void      CheckIfSinglet();

    //Event based functions, inline
    void BookkeepScatter(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, flavourpair out_flavours, double t);
    void BookkeepSplit(std::shared_ptr<Parton> part_in, flavourpair out_flavours, double t);
    void BookkeepMerge(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double t);

    void BookkeepOverdraw(double val, bool is_split);

    void ResetEventVariables();

  };
}

#endif
