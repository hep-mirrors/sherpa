#ifndef Full_Decay_Table_h
#define Full_Decay_Table_h

#include "Process_Group.H"
#include "Single_Process.H"
#include "Decay_Table.H"
#include <set>

namespace AMEGIC {
  class Full_Decay_Channel {
  private :
    ATOOLS::Decay_Channel   m_dec;
    Single_Process        * p_proc;
  public :
    Full_Decay_Channel(ATOOLS::Decay_Channel &);
    Full_Decay_Channel(const ATOOLS::Flavour);

    void AddDecayProduct(const ATOOLS::Flavour &);
    void Output();

    bool   CreateDecay();
    bool   InitProcess(Interaction_Model_Base *,Topology *);
    bool   CalculateWidth();
    void   SetWidth(double _w = -1.);

    double                  Width()           { return m_dec.Width(); }
    Process_Base          * GetProcessBase()  { return p_proc; }
    ATOOLS::Decay_Channel * GetDecayChannel() { return &m_dec; }
  };

  class Full_Decay_Table {
  private :
    bool                              m_isevaluated;
    double                            m_width;
    ATOOLS::Flavour                   m_flin;
    std::vector<Full_Decay_Channel *> m_channels;
    std::vector<Process_Group *>      m_decaymodes;
  public :
    Full_Decay_Table(const ATOOLS::Flavour);

    void   AddDecayChannel(ATOOLS::Decay_Channel *);
    void   AddDecayChannel(Full_Decay_Channel *);
    void   ArrangeDecays();
    bool   InitAllDecays(Interaction_Model_Base *,Topology *);

    void   CalculateWidths();
    void   Output();

    double                  Width()            { return m_width; }
    bool                    IsEvaluated()      { return m_isevaluated; }
    int                     NumberOfChannels() { return m_channels.size(); }
    ATOOLS::Decay_Channel * GetChannel(int);
  };
}


#endif
