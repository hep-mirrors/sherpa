#ifndef PHASIC_Main_Channel_Creator_H
#define PHASIC_Main_Channel_Creator_H

#include "ATOOLS/Org/CXXFLAGS.H"

namespace PHASIC {
  class Phase_Space_Handler;
  class Channel_Creator {
  private:
    Phase_Space_Handler * p_psh;
    
    bool CreateBeamIntegrator(const double & minalpha);
    bool CreateISRIntegrator(const double & minalpha);
    bool CreateFSRIntegrator(const double & minalpha);
  public:
    Channel_Creator(Phase_Space_Handler * psh);
    ~Channel_Creator();

    bool operator()();
  };
}
#endif
