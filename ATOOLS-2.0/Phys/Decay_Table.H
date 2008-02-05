#ifndef Decay_Table_h
#define Decay_Table_h

#include "Flavour.H"
#include "Message.H"
#include <set>

namespace ATOOLS {
  class Decay_Channel {
  protected :
    std::string m_processname;
    std::string m_metype, m_psfile;
    double      m_width, m_minmass;
    Flavour     m_flin;
    FlavourSet  m_flouts;
  public :
    Decay_Channel() {}
    Decay_Channel(const Flavour &);
    Decay_Channel(const Decay_Channel &);
    void AddDecayProduct(const Flavour &);
    void SetWidth(const double &);
    void Output() const;
    void SetProcessName(const std::string=std::string(""));
    void SetMEType(const std::string=std::string(""));
    void SetPSFile(const std::string=std::string(""));

    std::string ProcessName() const;
    std::string METype() const;
    std::string PSFile() const;
    double      Width() const;
    Flavour     GetDecaying() const;
    int         NumberOfDecayProducts() const;
    FlavourSet  GetDecayProducts() const;
    Flavour     GetDecayProduct(int) const;
    double      MinimalMass() const;
  };

  inline void Decay_Channel::AddDecayProduct(const Flavour & _flout) { 
    m_flouts.insert(_flout); 
    m_minmass += _flout.Mass();
  }
  inline void        Decay_Channel::SetWidth(const double & _width)           { m_width  = _width; }
  inline double      Decay_Channel::Width() const                             { return m_width;   }
  inline std::string Decay_Channel::ProcessName() const                       { return m_processname; }
  inline void        Decay_Channel::SetMEType(const std::string _metype)      { m_metype = _metype; }
  inline void        Decay_Channel::SetPSFile(const std::string _psfile)      { m_psfile = _psfile; }
  inline std::string Decay_Channel::METype() const                            { return m_metype; }
  inline std::string Decay_Channel::PSFile() const                            { return m_psfile; }
  inline Flavour     Decay_Channel::GetDecaying() const                       { return m_flin;    }
  inline int         Decay_Channel::NumberOfDecayProducts() const             { return m_flouts.size();  }
  inline double      Decay_Channel::MinimalMass() const                       { return m_minmass; }
  inline FlavourSet  Decay_Channel::GetDecayProducts() const                  { return m_flouts;  }
  inline Flavour     Decay_Channel::GetDecayProduct(int i) const {
    if (i<0 || i>(int)m_flouts.size()) {
      msg_Out()<<"WARNING in Decay_Channel::GetDecayProduct("<<i<<") :"<<std::endl
	       <<"   out of bounds, return no flavour."<<std::endl;
      return Flavour(kf::none);
    }
    for (FlSetConstIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) {
      if (i==0) return (*flit);
      --i;
    }
    return Flavour(kf::none);
  }


  class Decay_Table {
  protected :
    bool                          m_overwrite, m_smearing, m_fixdecay;
    std::string                   m_generator;
    double                        m_width;
    Flavour                       m_flin;
    FlavourSet                    m_flouts;
  private :
    std::vector <Decay_Channel *> m_channels;
    Decay_Channel *               p_selected;
    std::vector <Decay_Channel *> m_selectedchannels, m_selectedchannelsbar;
    std::vector <Decay_Channel *>::iterator m_seliter, m_seliterbar;
  public :
    Decay_Table(const Flavour);
    void AddDecayChannel(Decay_Channel *);
    void Reset();
    void Select(const int flag=0);
    void SetSelectedChannel(const FlavourSet &,const bool bar=false);
    void SetOverwrite()                      { m_overwrite = 1; }
    void SetSmearing()                       { m_smearing  = 1; }
    void SetWidthGenerator(std::string _gen) { m_generator = _gen; }
    void Output();
    void UpdateWidth();

    bool            FixedDecay()             { return m_fixdecay;  }
    bool            Overwrite()              { return m_overwrite; }
    bool            Smearing()               { return m_smearing;  }
    Flavour         Flav()                   { return m_flin;      }
    FlavourSet      GetDecayFlavs()          { return m_flouts;    }
    double          TotalWidth()             { return m_width;     }
    int             NumberOfDecayChannels()  { return m_channels.size(); } 
    double          Width(const int);
    double          Width(const FlavourSet);
    Decay_Channel * GetDecayChannel(const int);
    Decay_Channel * GetDecayChannel(const FlavourSet);
    Decay_Channel * GetOneDecayChannel();
    void            EraseDecayChannel(const int);
  };
}

#endif
