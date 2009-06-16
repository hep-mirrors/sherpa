#ifndef ATOOLS_Phys_Decay_Table_H
#define ATOOLS_Phys_Decay_Table_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include <set>
#include "ATOOLS/Phys/Decay_Channel.H"

namespace ATOOLS {
  class Decay_Table : public std::vector<Decay_Channel *> {
  protected :
    bool                          m_overwrite, m_smearing, m_fixdecay;
    double                        m_width;
    Flavour                       m_flin;
  public :
    Decay_Table(const Flavour);
    virtual ~Decay_Table();
    void AddDecayChannel(Decay_Channel *);
    void Output();
    void UpdateWidth();
    void ScaleToWidth();

    Flavour         Flav()                   { return m_flin;      }
    double          TotalWidth()             { return m_width;     }
    Decay_Channel * GetDecayChannel(const FlavourSet);
    void            EraseDecayChannel(const int);
    Decay_Channel * Select();

    friend std::ostream &operator<<(std::ostream &os, const Decay_Table &dt);
  };
}

#endif
