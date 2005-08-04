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
	std::string m_path, m_file;
	bool m_fulldecay, m_antidecays, m_createbooklet;
  public: 
    Decay_Table_Reader( std::string path=std::string("./"), 
						std::string file=std::string(""));
    int  FillDecayTable(ATOOLS::Decay_Table * dt=NULL);
    void FillInMatrixElementsAndPS( ATOOLS::Decay_Table *dt, 
		Channel_Map *chmap );
  }; // end of class Decay_Table_Reader
  
  /*!
	\file Decay_Table_Reader.H
	\brief Declares the class HADRONS::Decay_Table_Reader

	This code can be found in directory <tt>Main</tt>.
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
	Furthermore, this method reads decay channel file and stores its content in the private attribute.
  */	
  /*!
	\fn Decay_Table_Reader::FillDecayTable(ATOOLS::Decay_Table * dt=NULL)
	\brief Creates decay table with branching ratios and decay widths

	Uses data stored in the private attribute  to create a decay table, 
	i.e. decay channel, branching ratio, decay width, etc. The method
	returns the number of possible decay channels. The argument, finally, is a pointer on the complete decay table.
  */	
  /*!
	\fn Decay_Table_Reader::FillInMatrixElementsAndPS(ATOOLS::Decay_Table * dt=NULL)
	\brief Initialises phases spaces for each decay channel

	Reads decay table and initialises for every channel a phase space using 
	Hadron_Decay_Channel::InitialisePhaseSpace.
  */	


	 
	  
  
} // end of namespace HADRONS


#endif
