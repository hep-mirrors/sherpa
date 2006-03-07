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
#include "Spin_Density_Matrix.H"
#include "Spin_Correlation_Tensor.H"


using namespace ATOOLS;

namespace HADRONS {
  
  class Hadrons {

  private:
    std::string           m_path, m_file, m_constfile;
	Channel_Map         * p_channelmap;
	ATOOLS::Decay_Table * p_table;                 
	std::map<ATOOLS::kf::code,ATOOLS::Decay_Table *> 
	                    * p_decaymap;
    GeneralModel          m_md0;                // starter model
	
    void ReadInConstants();
    void ReadInDecayTables();
    ATOOLS::Decay_Table * InitialiseOneDecayTable(std::vector<std::string>);

	bool				  m_createbooklet;

  public: 
    Hadrons(
        std::string _path=std::string("./"), 
        std::string _file=std::string("HadronDecays.dat"),
        std::string _constfile=std::string("HadronConstants.dat") );
	
	bool                   FindDecay( const Flavour & );
	Hadron_Decay_Channel * ChooseDecayChannel();
	void                   ChooseDecayKinematics( Vec4D * _p, 
                                                  Hadron_Decay_Channel * _hdc,
                                                  ATOOLS::Spin_Density_Matrix * sigma );
	Spin_Density_Matrix    PerformDecay( Particle * part, 
		                                 Blob_List * blob_list, 
		                                 Particle_List * part_list,
                                         ATOOLS::Spin_Density_Matrix * sigma );

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
	  which stores information about the DT-files.
      Frequently used constants are written in \c HadronConstants.dat.
	- <b>DT-file, Decay-table-file</b>: each hadron has its own DT-file 
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

	\section sc The Implementation of Spin Correlations

    The usage with spin correlation from a programming point of view is a bit involved.
    A detailed explanation is given in the description of the class HADRONS::HD_ME_Base.

	
  */	
	
  /*!
	\file Hadrons.H
	\brief Declares the class HADRONS::Hadrons

	This code can be found in the directory <tt>Main</tt>.
  */	
  /*!
	\class Hadrons
	\brief Main class for hadron decay handling
  */	
  
  /*!
	\fn Hadrons::Hadrons(std::string _path=std::string("./"), std::string _file=std::string("HadronDecays.dat"), std::string _constfile=std::string("HadronConstants.dat") )
	\brief Constructor for the hadron decay handling

	This constructor calls Hadrons::ReadInConstants and Hadrons::ReadInDecayTables 
    to read frequently used constants \c HadronConstants.dat and the Hadron Decays file \c HadronDecays.dat.
  */
  /*!
	\fn Hadrons::ReadInDecayTables()
	\brief Reads decay tables and evaluates them
	
	It reads the Hadron decays file and extracts the filenames of the DT-files.
	These files however
	are evaluated in Hadrons::InitialiseOneDecayTable whose argument is a vector of strings which was created
	using the ATOOLS::Data_Reader::MatrixFromFile method.
  */	
  /*!
    \fn Hadrons::ReadInConstants()
    \brief Reads frequently used constants

    This methods reads a file like \c HadronConstants.dat (the actual name is stored in Hadrons::m_constfile
    in the Decaydata directory Hadrons::m_path). It reads the parameters and stores them in a string matrix
    of the form
    <TABLE>
    <TR> <TD>1st column</TD><TD>2nd column</TD><TD>3rd column</TD> </TR>
    <TR> <TD>name</TD><TD>=</TD><TD>value</TD> </TR>
    </TABLE>
    These will be the first setup for the model map HD_ME_Base::GeneralModel, which is used
    for passing parameters from the DC files into the ME class. See HD_ME_Base::GeneralModel for more details.
    */
  
  /*!
	\fn Hadrons::InitialiseOneDecayTable(std::vector<std::string>)
	\brief Evaluates a decay table of a hadron

	This Method uses the Decay_Table_Reader::FillDecayTable method to create the decay table of the considered
	hadron. It also calls Decay_Table_Reader::FillInMatrixElementsAndPS that actually reads the
    DT-file.
  */	
  /*!
    \fn Hadrons::FindDecay(const ATOOLS::Flavour & Decayer)
    
    Returns true if the Hadrons module can cope with the particle <var>Decayer</var>.
    */
  /*!
    \fn Hadrons::ChooseDecayChannel()
    \brief Chooses a decay channel according to branching ratios

    It chooses according to the given BR in a decay table a decay channel for the
    decaying hadron. Hereby, BR are normalised such that their sum equals to one (if necessary).
    The method returns a pointer on the chosen decay channel (Hadron_Decay_Channel class).

    For the K0 particle it chooses with a 50-50 probability either K0->KL or K0->KS.
    */
  /*!
    \fn Hadrons::ChooseDecayKinematics( ATOOLS::Vec4D * _p, Hadron_Decay_Channel * _hdc, ATOOLS::Spin_Density_Matrix * sigma ) 
    \brief Chooses the decay kinemtics

    This routine uses the trial-and-error method to obtain the kinematics of the decay products.
    That means, a phase space point is chosen and its value \f$|M|^2 \times dlips\f$ is compared with the maximum
    Hadron_Decay_Channel::m_max. According to this ratio the phase space point is rejected or not.
    In case of spin correlations the maximum needs to be overestimated (because the maximum was obtained
    during the uncorrelated integration HD_PS_Base::CalculateNormalisedWidth). 
    
    It calls the 
    Hadron_Decay_Channel::Differential( Vec4D * mom, Spin_Density_Matrix * sigma ) 
    method to obtain the value \f$|M|^2 \times dlips\f$. <var>mom</var> is a pointer on the momenta
    ans <var>sigma</var> a pointer on the spin density matrix. If there is none, <var>sigma</var> points
    on \c NULL.
    */
  /*!
    \fn Hadrons::PerformDecay( ATOOLS::Particle * part, ATOOLS::Blob_List  * blob_list, ATOOLS::Particle_List * part_list, ATOOLS::Spin_Density_Matrix * sigma )
    \brief Performs the decay of the particle <var>part</var>

    This is the main method for the decay of the particle <var>part</var>
    -# It choses a decay channel by calling Hadrons::ChooseDecayChannel.
    -# If this channel has the obtion \c NO_ANTI (Hadron_Decay_Channel::m_fulldecay) and the decaying particle is an anti-particle
       the methods does not continue.
    -# The kinematics of the daughters by calling Hadrons::ChooseDecayKinematics.
    -# A new blob is created.
    -# (if applicable) Builds up a spin correlation tensor using the amplitude tensor and index list 
       from the choses hadron decay channel (Hadron_Decay_Channel::p_ampls and Hadron_Decay_Channel::p_indices).
    -# Treat every daughter in a random order:
        - Hadrons::PerformDecay is called self-recursively
        - Spin correlation algorithm is fully implemented (see hep-ph/0110108)
        .
    .
    */

    
  
    

  


} // end of namespace


#endif
