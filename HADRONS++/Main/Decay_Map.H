#ifndef Decay_Map_H
#define Decay_Map_H

#include <string>
#include <vector>
#include "Flavour.H"
#include "Tools.H"

namespace HADRONS {
  class Hadron_Decay_Table;

  typedef std::pair<std::string,std::pair<double,double> >  BRPairString;
  typedef std::pair<ATOOLS::FlavourSet,std::pair<double,double> >   BRPairFlavourSet;
  typedef std::map<std::string,std::pair<double,double> >   BRMapString;
  typedef std::pair<double,double> DoublePair;
  typedef std::map<ATOOLS::Flavour,std::vector<Hadron_Decay_Table *> > FlDtVMap;
  typedef std::multimap<std::string, Hadron_Decay_Table*> SDtMMap;
  typedef std::multimap<std::string, Hadron_Decay_Table*>::iterator SDtMMapIt;

  class Decay_Map : public FlDtVMap {
  private:
    std::map<ATOOLS::Flavour,int> m_counters;
    GeneralModel m_startmd;
    std::string m_decaypath, m_decayfile, m_constfile;
    std::multimap<std::string, Hadron_Decay_Table*> m_fixed_tables;
    std::vector<Hadron_Decay_Table*> m_fixed_next_tables;

    std::vector<BRPairFlavourSet> GetInclusives(
        Hadron_Decay_Table * dt,
        BRMapString & brmap,
        ATOOLS::FlavourSet flset=ATOOLS::FlavourSet(),
        DoublePair br=DoublePair(1.,1.),
        bool eoi=true );
    void CreateBooklet();

  public: 
    Decay_Map(std::string path,
	      std::string file,
              std::string constfile );
    ~Decay_Map();
    void Read();
    void Initialise();
    void ReadFixedTables();
    void ReadInConstants();
    void SetHadronProperties();
    GeneralModel& StartModel() { return m_startmd; }
    Hadron_Decay_Table* FindDecay( const ATOOLS::Flavour & );
    void ResetCounters();
    void FixDecayTables(std::string table_id);
    void ClearFixedDecayTables();
  }; // end of class Decay_Map
  
  /*!
	\file Decay_Map.H
	\brief Declares the class HADRONS::Decay_Map

	This code can be found in the directory <tt>Main</tt>.
  */

  /*!
	\class Decay_Map
	\brief Tool to read decay tables
  */

  /*!
	\fn Decay_Map::Decay_Map( std::string path=std::string("./"), std::string file=std::string(""))
	\brief Constructor using ATOOLS::Data_Reader class
	
    Initializes the Decay_Map with the cut "|",
    the column separator ";" and the comment indicator "!". Additionally, "#" and "//" are interpreted as comments.
	Furthermore, this method reads decay channel file and stores its content in the private attribute 
    Decay_Map::m_helpsvv.
    
  */	
  /*!
	\fn Decay_Map::FillDecayTable(ATOOLS::Decay_Table * dt=NULL)
	\brief Creates decay table with branching ratios and decay widths

	Uses data stored in the private attribute  to create a decay table, 
	i.e. decay channel, branching ratio, decay width, etc. The method
	returns the number of possible decay channels. Finally, the argument is a pointer on the complete decay table.
  */	
  /*!
    \fn Decay_Map::ExtractFlavours(std::vector<int> &,std::string)
    \brief Converts an expression <a> {kf1,kf2,kf3,...} </a> into a vector with the kf-codes.
    */
  /*!
    \fn Decay_Map::FillInMatrixElementsAndPS( ATOOLS::Decay_Table *dt, Channel_Map *chmap, GeneralModel )
    \brief Creates pointer on Hadron_Decay_Channel

    This methods initialises a new hadron decay channel and puts a pointer on it. 
    It calls the Hadron_Decay_Channel::Hadron_Decay_Channel constructor for
    the initialisation of the decay channel and Hadron_Decay_Channel::InitialisePhaseSpace
    for the actual reading of the decay table.
    In case the decay table file needs to rewritten (that is the case when a DC file
    was generated automatically so that its name has to be inserted) it also 
    rewrites the decay table file with all known entries.
    */

  /*!
    \var Decay_Map::m_helpsvv
    
    A matrix of strings that contains a complete decay table with the following structure
    <TABLE>
    <TR> <TD>1st column</TD><TD>2nd column</TD><TD>3rd column</TD> </TR>
    <TR> <TD>decay products</TD><TD>branching ratio</TD><TD>name of DC file</TD> </TR>
    </TABLE>
    If a name of the DC file is \b not given, HADRONS creates a new DC file with a generic name.
    This DC file can be seen as template with the most basic entries (the intregrator is Rambo and 
    no parameters are given).
   */
  /*!
    \var Decay_Map::m_path

    Path to decay tables, e.g. \c TauDecays/, \c Rho0Decays/, etc.
   */
  /*!
    \var Decay_Map::m_file

    Name of the file with the decay table, e.g. \c TauDecays.dat, \c Rho0Decays.dat, etc.
    */



	 
	  
  
} // end of namespace HADRONS


#endif
