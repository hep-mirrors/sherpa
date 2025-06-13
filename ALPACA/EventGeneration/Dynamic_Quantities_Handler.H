#ifndef ALPACA_Event_Generation_Dynamic_Quantities_Handler_H
#define ALPACA_Event_Generation_Dynamic_Quantities_Handler_H

#include "ALPACA/Tools/Parton.H"
#include "ALPACA/Tools/HeavyIon_Parameters.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "MODEL/Main/Strong_Coupling.H"

#include <memory>
#include <vector>
#include <iostream>



namespace ALPACA {

  class Analysis_Handler;

  class Dynamic_Quantities_Handler {
  private :
    //General variables
    double                    m_alphaS;
    MODEL::Strong_Coupling *  p_alphaS;

    std::shared_ptr<std::list<std::shared_ptr<Parton>>> p_partons;
    std::shared_ptr<std::list<std::shared_ptr<Parton>>> p_removed_partons;

    std::shared_ptr<Analysis_Handler>                   p_analysis_handler;

    int                     m_timekeeper;
    bool                    m_only_gluons;

    //For mg2, mq2, Tstar
    int       m_N_inc;
    bool      m_iterative_m2_Tstar;
    double    m_m2_min_scale;

    //For PSD
    double    m_f_r;
    double    m_f_delta_p;
    bool      m_f_shell;
    int       m_f_N_max;

    double          Abs3(ATOOLS::Vec4D v);
    
  public:
    Dynamic_Quantities_Handler(std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_partons, 
                               std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_removed_partons);
    ~Dynamic_Quantities_Handler();

    void SetAnalysisHandler(std::shared_ptr<Analysis_Handler> ptr_analysis_handler);


    //Main functions
    double                    GetRunningAlphaS(double shat);
    double                    GetAlphaS();
    std::vector<double>       Getm2Tstar(std::shared_ptr<Parton> part1, std::shared_ptr<Parton> part2, double tau);
    double                    GetPSD(std::shared_ptr<Parton> part1, std::shared_ptr<Parton> part2, 
                                     double tau, ATOOLS::Vec4D p, bool isQuark, double Tstar);
    std::vector<std::pair<double, std::shared_ptr<Parton>>>  FindRecoilPartner(double tau, std::shared_ptr<Parton> part_1, 
                                                                               std::shared_ptr<Parton> part_2);

    //Misc. functions
    double          FindClosestR(ATOOLS::Vec4D x1, ATOOLS::Vec4D x2); 

    static bool sortM2(std::pair<double,std::pair<std::pair<ATOOLS::Vec4D,ATOOLS::Vec4D>,ATOOLS::Flavour>> i, std::pair<double,std::pair<std::pair<ATOOLS::Vec4D,ATOOLS::Vec4D>,ATOOLS::Flavour>> j) { return (i.first<j.first); }
    static bool sortDistList(std::pair<double,std::shared_ptr<Parton>> i, std::pair<double,std::shared_ptr<Parton>>j) { return (i.first<j.first); }
  };
}

#endif
