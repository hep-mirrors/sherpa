#ifndef USING__Hadrons
#define USING__Hadrons
#endif

#ifndef Hadrons_H
#define Hadrons_H

#include "CXXFLAGS.H"
#include "Decay_Table_Reader.H"
#include "Decay_Table.H"
#include <string>
#include <map>
#include "Vector.H"
#include "Hadron_Decay_Channel.H"
#include "Particle.H"
#include "Particle_List.H"
#include "Blob_List.H"


using namespace ATOOLS;

namespace HADRONS {
  
  class Hadrons {

  private:
    std::string           m_path, m_file;
	Channel_Map         * p_channelmap;
	ATOOLS::Decay_Table * p_table;                 
	std::map<ATOOLS::kf::code,ATOOLS::Decay_Table *> 
	                    * p_decaymap;
	
    void ReadInDecayTables();
    ATOOLS::Decay_Table * InitialiseOneDecayTable(std::vector<std::string>);

	bool				  m_createbooklet;

  public: 
    Hadrons(std::string _path=std::string("./"), std::string _file=std::string("HadronDecays.dat") );
	
	bool                   FindDecay( const Flavour & );
	Hadron_Decay_Channel * ChooseDecayChannel();
	void                   ChooseDecayKinematics( Vec4D * _p, Hadron_Decay_Channel * _hdc );
//	double                 Lifetime( double tau );
	void                   PerformDecay( Particle * part, 
		                                 Blob_List * blob_list, 
		                                 Particle_List * part_list );

	std::map<ATOOLS::kf::code,ATOOLS::Decay_Table *> * GetDecayMap() { return p_decaymap; }

	void CreateBookletNow();
  }; // end of class Hadrons
  
  /*!
	\mainpage

	The HADRONS package is the tool which deals with hadron decays. It provides a collection
	of ME and PS handlings. 

	\section files Dat-Files

	The idea is to save necessery information for the decays in files. Throughout this
	documentation three different terms are used each of which stands for a group of dat-files:
	- <b> Hadron decays file</b>: by default it is the file \c HadronDecays.dat, 
	  which stores information about the H-files.
	- <b>H-file, Hadron file</b>: each hadron has its own H-file 
	  which stores its decay modes and the
	  DC-filenames of the corresponding channel.
	- <b> DC-file, Decay channel file</b>: it has all necessery information about a complete
	  decay channel
	.

	\section libraries Libraries

	The HADRONS package has a couple of libraries, they are grouped by
	-# <b>PS Libraries</b>: tools to deal with the phase space of a certain decay.
	-# <b>ME Libraries</b>: tools to calculate the matrix element of a certain decay.
	.


	
  */	
	
  /*!
	\file Hadrons.H
	\brief Declares the class HADRONS::Hadrons

	This code can be found in directory <tt>Main</tt>.
  */	
  /*!
	\class Hadrons
	\brief Tool to deal with a decay table
  */	
  
  /*!
	\fn Hadrons::Hadrons(std::string _path=std::string("./"), std::string _file=std::string("HadronDecays.dat"))
	\brief Constructor for a hadron decay

	This constructor calls Hadrons::ReadInDecayTables to read the Hadron decays file.
  */
  /*!
	\fn Hadrons::ReadInDecayTables()
	\brief Reads decay tables and evaluates them
	
	It reads the Hadron decays file and extracts the filenames of the H-files.
	These files, however,
	are evaluated in Hadrons::InitialiseOneDecayTable whose argument is a vector of strings which was created
	using the ATOOLS::Data_Reader::MatrixFromFile method.
  */	
  /*!
	\fn Hadrons::InitialiseOneDecayTable(std::vector<std::string>)
	\brief Evaluates a decay table of a hadron

	This Method uses Decay_Table_Reader::FillDecayTable method to create the decay table of the considered
	hadron and has phase spaces for each channel initialised by calling Decay_Table_Reader::FillInMatrixElementsAndPS.
  */	

  


} // end of namespace


#endif
