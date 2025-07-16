#ifndef HADRONS_Main_Hadron_Decay_Map_H
#define HADRONS_Main_Hadron_Decay_Map_H

#include <string>
#include <vector>
#include "ATOOLS/Phys/Flavour.H"
#include "HADRONS++/Main/Tools.H"
#include "PHASIC++/Decays/Decay_Map.H"

namespace HADRONS {
  class Hadron_Decay_Table;
  class Mixing_Handler;

  typedef std::multimap<std::string, Hadron_Decay_Table*> SDtMMap;
  typedef std::multimap<std::string, Hadron_Decay_Table*>::iterator SDtMMapIt;

  /**
   * Building on the general class Decay_Map this class implements the features
   * specific for hadron decay maps, e.g. reading decay tables from files.
   * Furthermore it manages the global constants defined in HadronConstants.dat.
   */
  class Hadron_Decay_Map : public PHASIC::Decay_Map {
  private:

    /**
     * A model that is filled from HadronConstants.dat and used as starting
     * point for all decay channel specific models.
     */
    GeneralModel m_startmd;

    std::multimap<std::string, Hadron_Decay_Table*> m_fixed_tables;
    std::vector<Hadron_Decay_Table*> m_fixed_next_tables;

    Mixing_Handler      * p_mixinghandler;

  public: 
    /** 
     * @param ms The mass selector
     */
    Hadron_Decay_Map(const ATOOLS::Mass_Selector* ms);
    ~Hadron_Decay_Map();
    /**
     * Reads all decay tables and create their decay channels and initialises them.
     */
    void Read(ATOOLS::Scoped_Settings& s);
    /**
     * Reads fixed decay tables from FixedDecays.dat which can later be
     * programmatically set to be used for the next decay of a certain decayer,
     * typically some kind of signal decay
     * (useful within frameworks like LHCb)
     */
    void ReadFixedTables();
    /**
     * Reads partonic decay data.
     */
    void ReadInPartonicDecays(const ATOOLS::Flavour & decflav, ATOOLS::Scoped_Settings& s);
    /** 
     * 
     * @return The global hadron decay parameter model.
     */
    inline GeneralModel& StartModel() { return m_startmd; }
    /**
     * Creates a Latex booklet containing all decay tables and information
     * about available matrix elements and currents, in hadrons.tex
     */
    void CreateBooklet();
    /**
     * Activate a previously read fixed decay table 
     */
    void FixDecayTables(std::string table_id);
    /**
     * Deactivate the previously fixed decay table
     */
    void ClearFixedDecayTables();

    virtual PHASIC::Decay_Table* FindDecay(const ATOOLS::Flavour& decayer);

    //@{
    /// Getter/setter method.
    inline const Mixing_Handler* MixingHandler() const {return p_mixinghandler;}
    inline void SetMixingHandler(Mixing_Handler* mix) { p_mixinghandler = mix; }
    //@}
  };

  /**
     \package HADRONS
     The HADRONS package is Sherpa's native hadron and tau decay handler.

     \section decaydata Decaydata files

     All decaytables and information about decay channels is saved in plain
     text files located in SHERPA/Run/Decaydata.

     Three different levels of dat-files are used:
     - <b>HadronDecays.dat</b>: Stores a list of all available decay tables.
     (represented by the Hadron_Decay_Map class).
     - <b>Decay table files</b>: Each hadron has its own directory including
     a file called Decays.dat which is represented by the
     Hadron_Decay_Table class.
     - <b>Decay channel files</b>: In the same directory you can find individual
     files for each decay channel (represented by the Hadron_Decay_Channel
     class).
  */

}


#endif
