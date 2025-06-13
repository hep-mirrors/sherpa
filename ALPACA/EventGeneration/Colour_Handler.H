#ifndef ALPACA_Event_Generation_Colour_Handler_H
#define ALPACA_Event_Generation_Colour_Handler_H

#include "ALPACA/Tools/Parton.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Phys/Flavour.H"

#include <memory>
#include <vector>
#include <map>
#include <list>

namespace ALPACA {

  class Colour_Handler {
  private:
    std::shared_ptr<std::list<std::shared_ptr<Parton>>> p_partons;
    std::shared_ptr<std::list<std::pair<int,int>>>      p_flow_backtrack; //List of colours to keep track of where the colour partner is in a blob not converted and sent to ALPACA
    
  public:
    Colour_Handler(std::shared_ptr<std::list<std::shared_ptr<Parton>>> ptr_partons, std::shared_ptr<std::list<std::pair<int,int>>> ptr_flow_backtrack);
    ~Colour_Handler();

    //General colour update functions for scattering/splitting/merging
    void UpdateColoursScatter(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                              ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                              double t, double mg2, double mf2);
    
    void UpdateColoursSplit(std::shared_ptr<Parton> part_in,
                            std::shared_ptr<Parton> part_out_1, std::shared_ptr<Parton> part_out_2);

    void UpdateColoursMerge(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                            std::shared_ptr<Parton> part_out);

    //Function to backtrack and change colour indices
    void BackTrackColour(int flow_id_1, int flow_id_2, 
                         std::shared_ptr<Parton> part_exc_1, std::shared_ptr<Parton> part_exc_2);
    
    //Infinitely soft processes to switch colour indices, to avoid ending up with colour neutral gluons in mergings
    void InfSoftGGtoG(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2);

    void InfSoftQQbtoG(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2);

    //Process specific colour updates for scatterings
    bool UpdateColours_qiqi_qiqi(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                 ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                                 double t, double mg2);
    bool UpdateColours_qiqj_qiqj(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                 ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2);
    bool UpdateColours_qiqbi_qiqbi(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                   ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                                   double t, double mg2);
    bool UpdateColours_qiqbi_qjqbj(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                   ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2);
    bool UpdateColours_qiqbi_gg(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                                double t, double mg2, double mf2);
    bool UpdateColours_qig_qig(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                               ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                               double t, double mg2, double mf2);
    bool UpdateColours_gg_gg(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                             ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                             double t, double mg2);
    bool UpdateColours_gg_qiqbi(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, 
                                ATOOLS::Flavour flav_out_1, ATOOLS::Flavour flav_out_2, 
                                double t, double mg2, double mf2);
  };
}
#endif
