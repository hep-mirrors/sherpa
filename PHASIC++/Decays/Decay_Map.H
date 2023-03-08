#ifndef Decay_Map_h
#define Decay_Map_h

#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Channel.H"

namespace ATOOLS {
  class Mass_Selector;
}

namespace PHASIC {

  class FlavourComp {
    const ATOOLS::Mass_Selector* p_ms;
  public:
    FlavourComp(const ATOOLS::Mass_Selector* ms) : p_ms(ms) {}
    bool operator()(const ATOOLS::Flavour& f1, const ATOOLS::Flavour& f2) const;
  };

  class Decay_Map : 
    public std::map<ATOOLS::Flavour, Decay_Table *, FlavourComp> {
  protected:
    /**
     * The Mass_Selector knows which type of mass is to be used for all decays
     * associated with this decay map. This can be the HadMass() if the decays
     * belong to the non-perturbative regime of the event, or the Mass() if
     * they belong to the perturbative part. 
     */
    const ATOOLS::Mass_Selector* p_ms;

  public :
    /**
     * @param ms Mass_Selector, cf. p_ms
     * 
     * Constructor
     */
    Decay_Map(const ATOOLS::Mass_Selector* ms);
    virtual ~Decay_Map();

    /**
     * @param decayer The decaying particle flavour.
     * 
     * @return Whether this decay map contains decay table(s) for this flavour
     */
    bool Knows(const ATOOLS::Flavour& decayer);

    /**
     * @param decayer The decaying particle flavour.
     * 
     * @return The currently active decay table.
     *
     * This takes into account that multiple tables can exist for a Flavour
     * of which the n'th will be chosen for the n'th particle of that Flavour
     * in this event. (cf. Decay_Map::m_counters)
     */
    virtual Decay_Table* FindDecay(const ATOOLS::Flavour& decayer);

    /**
     * @param name Name of the wanted decay channel in the format a->b c ...
     * 
     * @return The decay channel and its decay table.
     */
    //std::pair<Decay_Table*, Decay_Channel*> FindDecayChannel(std::string name);

    /**
     * Resets the counters (cf. Decay_Map::m_counters) for flavour types.
     * Typically used at the end of one event.
     */
    void ResetCounters();

    friend std::ostream &operator<<(std::ostream &os, const Decay_Map &dm);
  };
}

#endif
