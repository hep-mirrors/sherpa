#ifndef Decay_Table_Reader_H
#define Decay_Table_Reader_H

#include "Decay_Table.H"
#include <string>
#include <vector>
#include "Hadron_Decay_Channel.H"

typedef std::map<ATOOLS::Decay_Channel *, HADRONS::Hadron_Decay_Channel *> Channel_Map;

namespace HADRONS {
  
  class Decay_Table_Reader {
  private:
    std::vector<std::vector<std::string> > m_helpsvv;
    bool ExtractFlavours(std::vector<int> &,std::string);
    void ExtractBRInfo( std::string entry, std::string & sbr, std::string & sdbr, std::string & sorigin );
	std::string m_path, m_file;
	bool m_fulldecay, m_antidecays, m_createbooklet;
  public: 
    Decay_Table_Reader( std::string path=std::string("./"), 
						std::string file=std::string(""));
    int  FillDecayTable(ATOOLS::Decay_Table * dt=NULL);
    void FillInMatrixElementsAndPS( ATOOLS::Decay_Table *dt, 
		Channel_Map *chmap, GeneralModel );
  }; // end of class Decay_Table_Reader
  
  /*!
	\file Decay_Table_Reader.H
	\brief Declares the class HADRONS::Decay_Table_Reader

	This code can be found in the directory <tt>Main</tt>.
  */

  /*!
	\class Decay_Table_Reader
	\brief Tool to read decay tables
  */

  /*!
	\fn Decay_Table_Reader::Decay_Table_Reader( std::string path=std::string("./"), std::string file=std::string(""))
	\brief Constructor using ATOOLS::Data_Reader class
	
    Initializes the Decay_Table_Reader with the cut "|", 
    the column separator ";" and the comment indicator "!". Additionally, "#" and "//" are interpreted as comments.
	Furthermore, this method reads decay channel file and stores its content in the private attribute 
    Decay_Table_Reader::m_helpsvv.
    
  */	
  /*!
	\fn Decay_Table_Reader::FillDecayTable(ATOOLS::Decay_Table * dt=NULL)
	\brief Creates decay table with branching ratios and decay widths

	Uses data stored in the private attribute  to create a decay table, 
	i.e. decay channel, branching ratio, decay width, etc. The method
	returns the number of possible decay channels. Finally, the argument is a pointer on the complete decay table.
  */	
  /*!
    \fn Decay_Table_Reader::ExtractFlavours(std::vector<int> &,std::string)
    \brief Converts an expression <a> {kf1,kf2,kf3,...} </a> into a vector with the kf-codes.
    */
  /*!
    \fn Decay_Table_Reader::FillInMatrixElementsAndPS( ATOOLS::Decay_Table *dt, Channel_Map *chmap, GeneralModel )
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
    \var Decay_Table_Reader::m_helpsvv
    
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
    \var Decay_Table_Reader::m_path

    Path to decay tables, e.g. \c TauDecays/, \c Rho0Decays/, etc.
   */
  /*!
    \var Decay_Table_Reader::m_file

    Name of the file with the decay table, e.g. \c TauDecays.dat, \c Rho0Decays.dat, etc.
    */



	 
	  
  
} // end of namespace HADRONS


#endif
