#ifndef Hadron_Wave_Function_H
#define Hadron_Wave_Function_H

#include "Flavour.H"
#include <map>
#include <string>


namespace AHADIC {
  typedef std::pair<ATOOLS::Flavour,ATOOLS::Flavour> FlavPair;
  typedef std::map<FlavPair *,double>                WFcomponent;
  typedef WFcomponent::iterator                      WFcompiter;


  class Hadron_Wave_Function {
  private:
    ATOOLS::Flavour m_hadron;
    int             m_spin, m_kfcode;
    WFcomponent     m_waves;
    bool            m_barrable;
  public:
    Hadron_Wave_Function();
    Hadron_Wave_Function(const ATOOLS::Flavour &);    
    ~Hadron_Wave_Function();

    void AddToWaves(FlavPair *,double);
    Hadron_Wave_Function * Anti();
    double Weight(ATOOLS::Flavour,ATOOLS::Flavour);

    ATOOLS::Flavour SetFlavour(ATOOLS::Flavour & flav) { m_hadron = flav;    }
    ATOOLS::Flavour GetFlavour() const                 { return m_hadron;    }
    int  KfCode() const                                { return m_kfcode;    }
    void SetKfCode(const int _kfcode)                  { m_kfcode = _kfcode; }
    int  Spin() const                                  { return m_spin;      }
    void SetSpin(const int _spin)                      { m_spin   = _spin;   }
    WFcomponent * GetWaves()                           { return &m_waves;    }
    bool ExistAnti()                                   { return m_barrable;  }

    friend std::ostream & operator<<(std::ostream & s, Hadron_Wave_Function & wf);
  };

  typedef std::map<ATOOLS::Flavour,Hadron_Wave_Function *> Hadron_WF_Map; 
  typedef Hadron_WF_Map::iterator                          Hadron_WF_Miter;
}

#endif
