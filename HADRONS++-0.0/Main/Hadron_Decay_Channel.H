#ifndef Hadron_Decay_Channel_H
#define Hadron_Decay_Channel_H

#include "Integrable_Base.H"
#include "HD_ME_Base.H"
#include "HD_PS_Base.H"
#include "Decay_Table.H"
#include <string>

namespace HADRONS {
  class Hadron_Decay_Channel : public PHASIC::Integrable_Base {
  private:
    ATOOLS::Decay_Channel * p_dc;
    HD_ME_Base            * p_me;
    HD_PS_Base            * p_ps;
    std::string             m_metype;
    void                    CalculateNormalisedWidth();
  public: 
    Hadron_Decay_Channel(ATOOLS::Decay_Channel * _dc);
    ~Hadron_Decay_Channel();
    void   InitialisePhaseSpace(std::vector<std::string> &);
    double Differential();
    double Differential(const ATOOLS::Vec4D *momenta) { return 0.; }
    double Differential2()                            { return 0.; }

    void   SwapInOrder();
    void   RestoreInOrder();
    void   DeSelect();
    bool   SelectOne();
    bool   ReSelect(int);
    bool   SelectOneFromList();
    HD_ME_Base * GetME() { return p_me; }
    HD_PS_Base * GetPS() { return p_ps; }
  };
}


#endif
