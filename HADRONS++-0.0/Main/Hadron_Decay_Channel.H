#ifndef Hadron_Decay_Channel_H
#define Hadron_Decay_Channel_H

#include "Integrable_Base.H"
#include "HD_ME_Base.H"
#include "HD_PS_Base.H"
#include "Decay_Table.H"
#include <string>
#include "Vector.H"

#define HBAR_C 197.326968e-15	// in [GeV mm] used for decay position

namespace HADRONS {
  class Hadron_Decay_Channel : public PHASIC::Integrable_Base {
  private:
    ATOOLS::Decay_Channel * p_dc;
    HD_ME_Base            * p_me;
    HD_PS_Base            * p_ps;
    std::string             m_metype;
    void                    CalculateNormalisedWidth();
	std::string				m_path;
	short					m_fulldecay;	
	bool					m_createbooklet;
	std::string				m_channelname, m_chnamenumbers;
  public: 
    Hadron_Decay_Channel( ATOOLS::Decay_Channel * _dc, std::string _path );
    ~Hadron_Decay_Channel();
    bool   InitialisePhaseSpace(std::vector<std::string> & PStype, GeneralModel startmd );
    void   WriteModelOnScreen( GeneralModel _locmd );
    void   SwapInOrder()    {}
    void   RestoreInOrder() {}
    void   DeSelect()       {}
    bool   SelectOne()      { return 1; }
	bool   ReSelect( int )  { return 1; }
	bool   SelectOneFromList() { return 1; }
    double Differential();
    double Differential(double&);
    double Differential(const ATOOLS::Vec4D *momenta) { return 0.; }
    double Differential2()                            { return 0.; }

	void   SetFullDecay( short a ) { m_fulldecay = a; }
	short  FullDecay() { return m_fulldecay; }
	void   SetCreateBooklet( bool a ) { m_createbooklet = a; }
	bool   CreateBooklet() { return m_createbooklet; }

    HD_ME_Base            * GetME() { return p_me; }
    HD_PS_Base            * GetPS() { return p_ps; }
	ATOOLS::Vec4D         * Momenta() { return p_momenta; }		// not as const
	ATOOLS::Vec4D           Momentum( int i ) { return p_momenta[i]; }
	ATOOLS::Decay_Channel * DecayChannel() { return p_dc; }

	std::string				ChannelName()     { return m_channelname; }
	std::string				ChannelNameNumbers() { return m_chnamenumbers; }
	double                  MeanLifetime()    { return HBAR_C/p_ps->Result(); } // in [mm]
	double                  InMass()          { return (p_dc->GetDecaying()).Mass(); }
	
  }; // end of class

  /*!
	\file Hadron_Decay_Channel.H
	\brief Declares the class HADRONS::Hadron_Decay_Channel

	This code can be found in directory \c Main.
  */	

  /*!
	\class Hadron_Decay_Channel
	\brief Tool to manage a certain hadronic decay channel

	This class is a subclass of PHASIC::Integrable_Base.
  */	
  
  /*!
	\fn Hadron_Decay_Channel::Hadron_Decay_Channel( ATOOLS::Decay_Channel * _dc )
	\brief Constructor for a hadronic decay channel

	Sets initial values such as result path and file name, decay name, flavours of involved particles
	and uses HD_ME_Selector::GetME to map the corresponding decay's ME.
  */
  /*!
	\fn Hadron_Decay_Channel::InitialisePhaseSpace(std::vector<std::string> & PStype )
	\brief Method to have the phase space initialised.

	Uses HD_PS_Base::HD_PS_Base in order to set up the phase space correctly and calls 
	HD_PS_Base::Initialise for initialisation if the PS generator has to be initialised,
	i.e. no settings from a last run is written in file. If there are, these integration channel settings
	are used for PS integration.
  */	
  /*!
	\fn Hadron_Decay_Channel::Differential()
	\brief Returns differential rate \f$d\Gamma \times dlips\f$.

	Uses HS_PS_Base methods to get a connection to the choosen phase space generator
	(such as multichannel, Rambo, etc.) and make it generate a point and calculate the 
	weight factor. Finally, HD_ME_Base::HD_ME_Base returns its value of the ME so that the return value is obtained.
	The particle decays from rest.
  */
  /*!
	\var Hadron_Decay_Channel::p_me
	Pointer on ME.
  */	
  /*!
	\var Hadron_Decay_Channel::p_ps
	Pointer on PS generator.
  */
  /*!
	\var Hadron_Decay_Channel::p_dc
	Pointer on decay channel.
  */ 
  /*!
	\var Hadron_Decay_Channel::m_metype
	Type of ME gained from ATOOLS::Decay_Channel class.
  */ 
  /*!
    /var Hadron_Decay_Channel::m_fulldecay
	Has information about decay: syntax "(binary number) = decimal number"
	- (00) = 0 ... no full decay, anti particle does not decay
	- (01) = 1 ... no full decay, anti particle does decay
	- (10) = 2 ... full decay, anti particle does not decay
	- (11) = 3 ... full decay, anti particle does decay
	.
  */	
	

  
  
} // end of namespace


#endif
