#ifndef Decay_Table_h
#define Decay_Table_h

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include <set>

namespace ATOOLS {
  class Blob;
  class Mass_Selector;
}

namespace PHASIC {
  class Decay_Channel;

  class Decay_Table : public std::vector<Decay_Channel *> {
  protected :
    double                        m_activewidth, m_totalwidth;
    ATOOLS::Flavour               m_flin;
    /**
     * The Mass_Selector knows which type of mass is to be used for all decays
     * associated with this decay table. This can be the HadMass() if the decays
     * belong to the non-perturbative regime of the event, or the Mass() if
     * they belong to the perturbative part. 
     */
    const ATOOLS::Mass_Selector* p_ms;


  public :
    Decay_Table(const ATOOLS::Flavour, const ATOOLS::Mass_Selector* ms);
    virtual ~Decay_Table();
    void AddDecayChannel(Decay_Channel *);
    void RemoveDecayChannel(size_t i);
    void Output();
    void UpdateWidth();
    void UpdateWidth(Decay_Channel * dc,const double &width);

    const ATOOLS::Flavour& Flav() const      { return m_flin;      }
    const double &  TotalWidth() const       { return m_totalwidth; }
    const double &  ActiveWidth() const      { return m_activewidth; }
    Decay_Channel * GetDecayChannel(const ATOOLS::Flavour_Vector&) const;
    void            EraseDecayChannel(const int);
    Decay_Channel * Select() const;

    friend std::ostream &operator<<(std::ostream &os, const Decay_Table &dt);
  };
}

#endif
