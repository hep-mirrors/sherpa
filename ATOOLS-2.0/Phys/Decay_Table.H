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
    double     m_width;
    Flavour    m_flin;
    FlavourSet m_flouts;
  public :
    Decay_Channel() {}
    Decay_Channel(const Flavour &);
    Decay_Channel(const Decay_Channel &);
    void AddDecayProduct(const Flavour &);
    void SetWidth(const double &);
    void Output();

    double     Width();
    Flavour    GetDecaying();
    FlavourSet GetDecayProducts();
  };

  inline void Decay_Channel::AddDecayProduct(const Flavour & _flout) { m_flouts.insert(_flout); }
  inline void Decay_Channel::SetWidth(const double & _width)       { m_width = _width; }
  inline double Decay_Channel::Width()                             { return m_width;   }
  inline Flavour Decay_Channel::GetDecaying()                      { return m_flin;    }
  inline FlavourSet Decay_Channel::GetDecayProducts()              { return m_flouts;  }


  class Decay_Table {
  protected :
    bool                          m_overwrite;
    std::string                   m_generator;
    double                        m_width;
    Flavour                       m_flin;
  private :
    std::vector <Decay_Channel *> m_channels;
  public :
    Decay_Table(const Flavour);
    void AddDecayChannel(Decay_Channel *);
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
  };
}

#endif
