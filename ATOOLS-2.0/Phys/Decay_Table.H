#ifndef Decay_Table_h
#define Decay_Table_h

#include "Flavour.H"
#include <set>

namespace ATOOLS {
  class Flavour_Sorting_Criterion {
  public :
    bool operator() (const Flavour & _fl1,const Flavour & _fl2) {
      int fl1 = _fl1.Kfcode(), fl2 = _fl2.Kfcode();
      if (fl1<fl2) return 0;
      if (fl1==fl2 && _fl1.IsAnti() && !(_fl2.IsAnti())) return 0;
      return 1;
    }
  };

  typedef std::set<Flavour,Flavour_Sorting_Criterion> FlavourSet;
  typedef FlavourSet::iterator FlSetIter;

  class Decay_Channel {
  protected :
    std::string m_processname;
    double      m_width, m_minmass;
    Flavour     m_flin;
    FlavourSet  m_flouts;
  public :
    Decay_Channel() {}
    Decay_Channel(const Flavour &);
    Decay_Channel(const Decay_Channel &);
    void AddDecayProduct(const Flavour &);
    void SetWidth(const double &);
    void Output();
    void SetProcessName(const std::string &);

    std::string ProcessName();
    double      Width();
    Flavour     GetDecaying();
    int         NumberOfDecayProducts();
    FlavourSet  GetDecayProducts();
    Flavour     GetDecayProduct(int);
    double      MinimalMass();
  };

  inline void Decay_Channel::AddDecayProduct(const Flavour & _flout) { 
    m_flouts.insert(_flout); 
    m_minmass += _flout.Mass();
  }
  inline void        Decay_Channel::SetWidth(const double & _width)           { m_width = _width; }
  inline double      Decay_Channel::Width()                                   { return m_width;   }
  inline void        Decay_Channel::SetProcessName(const std::string & _name) { m_processname = _name; }
  inline std::string Decay_Channel::ProcessName()                             { return m_processname; }
  inline Flavour     Decay_Channel::GetDecaying()                             { return m_flin;    }
  inline int         Decay_Channel::NumberOfDecayProducts()                   { return m_flouts.size();  }
  inline double      Decay_Channel::MinimalMass()                             { return m_minmass; }
  inline FlavourSet  Decay_Channel::GetDecayProducts()                        { return m_flouts;  }
  inline Flavour     Decay_Channel::GetDecayProduct(int i) {
    if (i<0 || i>(int)m_flouts.size()) {
      std::cout<<"Error in Decay_Channel::GetDecayProduct("<<i<<") :"<<std::endl
	       <<"   out of bounds, return no flavour."<<std::endl;
      return Flavour(kf::none);
    }
    for (FlSetIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) {
      if (i==0) return (*flit);
      --i;
    }
    return Flavour(kf::none);
  }


  class Decay_Table {
  protected :
    bool                          m_overwrite;
    std::string                   m_generator;
    double                        m_width;
    Flavour                       m_flin;
  private :
    std::vector <Decay_Channel *> m_channels;
    Decay_Channel *               p_selected;
  public :
    Decay_Table(const Flavour);
    void AddDecayChannel(Decay_Channel *);
    void Select();
    void Output();
    void SetOverwrite()                      { m_overwrite = 1; }
    void SetWidthGenerator(std::string _gen) { m_generator = _gen; }

    bool            Overwrite()              { return m_overwrite; }
    Flavour         Flav()                   { return m_flin; }
    double          TotalWidth()             { return m_width; }
    int             NumberOfDecayChannels()  { return m_channels.size(); } 
    double          Width(const int);
    double          Width(const FlavourSet);
    Decay_Channel * GetDecayChannel(const int);
    Decay_Channel * GetDecayChannel(const FlavourSet);
    Decay_Channel * GetOneDecayChannel()     { return p_selected; }
  };
}

#endif
