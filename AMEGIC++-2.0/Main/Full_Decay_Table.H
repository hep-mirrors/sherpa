#ifndef Full_Decay_Table_h
#define Full_Decay_Table_h

#include "Process_Group.H"
#include "Single_Process.H"
#include "Decay_Table.H"
#include <set>

namespace AMEGIC {
  class Full_Decay_Channel {
  private :
    ATOOLS::Decay_Channel * p_dec;
    Single_Process        * p_proc;
  public :
    Full_Decay_Channel(ATOOLS::Decay_Channel *);
    Full_Decay_Channel(const ATOOLS::Flavour);

    void AddDecayProduct(const ATOOLS::Flavour &);
    void Output();

    bool   CreateDecay();
    bool   InitProcess(MODEL::Model_Base *,Topology *);
    bool   CalculateWidth();
    bool   OneEvent(double);
    void   SetWidth(double _w = -1.);

    double                  Width()           { return p_dec->Width(); }
    Process_Base          * GetProcessBase()  { return p_proc; }
    ATOOLS::Decay_Channel * GetDecayChannel() { return p_dec; }
    int                     NIn();
    int                     NOut();
    std::string             Name();
    const ATOOLS::Flavour * Flavours();            
    const ATOOLS::Vec4D   * Momenta();          
    int                     NumberOfDiagrams();
    Point                 * Diagram(int i);
  };
  inline bool              Full_Decay_Channel::OneEvent(double _mass) { return p_proc->OneEvent(_mass); }
  inline int               Full_Decay_Channel::NIn()                  { return p_proc->NIn(); }           
  inline int               Full_Decay_Channel::NOut()                 { return p_proc->NOut(); }          
  inline std::string       Full_Decay_Channel::Name()                 { return p_proc->Name(); }          
  inline const ATOOLS::Flavour * Full_Decay_Channel::Flavours()       { return p_proc->Flavours(); }
  inline const ATOOLS::Vec4D   * Full_Decay_Channel::Momenta() { 
    if (p_proc->Partner()==p_proc) return p_proc->Momenta();
    return p_proc->Partner()->Momenta();
  }          
  inline int               Full_Decay_Channel::NumberOfDiagrams()     { return p_proc->NumberOfDiagrams(); }
  inline Point *           Full_Decay_Channel::Diagram(int i)         { return p_proc->Diagram(i); }


  class Full_Decay_Table {
  private :
    double                            m_width;
    ATOOLS::Flavour                   m_flin;
    bool                              m_isevaluated,m_overwrite;
    std::vector<Full_Decay_Channel *> m_channels;
    std::vector<Process_Group *>      m_decaymodes;
  public :
    Full_Decay_Table(const ATOOLS::Flavour,bool=true);

    void   AddDecayChannel(ATOOLS::Decay_Channel *);
    void   AddDecayChannel(Full_Decay_Channel *);
    void   ArrangeDecays();
    bool   InitAllDecays(MODEL::Model_Base *,Topology *);

    void   CalculateWidths();
    void   Output();

    double                  Width()            { return m_width; }
    bool                    IsEvaluated()      { return m_isevaluated; }
    int                     NumberOfChannels() { return m_channels.size(); }
    ATOOLS::Decay_Channel * GetChannel(int);
    Full_Decay_Channel    * GetFullChannel(int);
  };
}


#endif
