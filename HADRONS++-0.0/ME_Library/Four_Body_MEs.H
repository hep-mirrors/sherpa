#ifndef Four_Body_MEs_H
#define Four_Body_MEs_H

#include "Hadron_Decay_ME.H"

namespace HADRONS {
  class P_2PFF : public Hadron_Decay_ME {
  public:
    P_2PFF(int,int,ATOOLS::Flavour *,
	   std::string _met=std::string("P_2PFF"));
    double operator()(const ATOOLS::Vec4D *);
  };
}


#endif
