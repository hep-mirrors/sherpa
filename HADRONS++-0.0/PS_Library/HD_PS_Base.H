#ifndef HD_PS_Base_H
#define HD_PS_Base_H

#include "Single_Channel.H"
#include "Multi_Channel.H"
#include "HD_ME_Base.H"  
#include "Flavour.H"
#include "Vector.H"
#include <string>

namespace HADRONS {


  class Hadron_Decay_Channel;

  class HD_Channel_Selector {
  public:
    PHASIC::Single_Channel * GetChannel(
		int nin, int nout, const ATOOLS::Flavour * flavs, std::string name);
  }; // end of class
  
  class HD_PS_Base : public PHASIC::Multi_Channel {
	void Initialize( GeneralModel& );
  private:
    HD_Channel_Selector  * p_channelselector;
    Hadron_Decay_Channel * p_hdc;
    std::string            m_path, m_file;
    double                 m_res,m_error,m_max,m_flux;
	bool                   m_read_dcfile;

    bool     Construct( GeneralModel & );
    void     AddChannel(std::string,double=1.);
    void     CalculateNormalisedWidth();
    bool     ReadIn();
    bool     WriteOut();
  public:
    HD_PS_Base(
		Hadron_Decay_Channel * hdc, 
		std::string _path,
		std::vector<std::string> & _pst, 
		bool & mustinit, 
		struct GeneralModel &,
		bool );
    virtual ~HD_PS_Base();
    void     Initialise();
	inline double Maximum() { return m_max; }
  }; // end of class

  /*!
	\file HD_PS_Base.H
	\brief Declares the classes HADRONS::HD_Channel_Selector and HADRONS::HD_PS_Base

	This file can be found in directory \c PS_Library.
  */

  /*!
	\class HD_Channel_Selector
	\brief Tool to select an integration channel
  */
  /*!
	\fn HD_Channel_Selector::GetChannel(int nin, int nout, const ATOOLS::Flavour * flavs, std::string name)
	\brief Returns pointer on a PS generator

	#in must be 1 and the #out greater equal than 2. The method's argument <var>name</var> 
	carries the information what kind
	of decay is to be considered (e.g. isotropic, Dalitz, etc.). After some sanity checks with the number of decay
	products the method returns the pointer on the corresponding PS generator such as Rambo, Iso2Channel, or Dalitz.
  */	


  /*!
	\class HD_PS_Base
	\brief Tool to set up the connection to PS generator 

	This class is a subclass of PHASIC::MultiChannel
  */	
  /*!
	\fn HD_PS_Base::HD_PS_Base(Hadron_Decay_Channel * hdc, std::vector<std::string> & _pst, bool & mustinit)
	\brief Constructor and initialisatio of private attributes

	The constructor calls HD_PS_Base::Construct to set up the PS. Its return value is referenced to the
	argument <var>mustinit</var>.
	Afterwards, the internal (private) pointer
	p_channelselector is deleted.
  */	
  /*!
	\fn HD_PS_Base::Construct
	\brief Constructs the PS according to given information

	This method reads the DC file, such as \f$\eta \to 2\gamma\f$, and
	does now extract stored information about
	  - Channels: kind of decay (isotropic, ...) and a weight factor (if given, if not: it's set to 1). 
	    These are piped 
	    into HD_PS_Base::AddChannel which is tasked with to set pointers right.
	  - Dalitz-Paramters: 5 Dalitz parameters which are processed in HD_ME_Base::SetDalitzParameters.
	  - Result: result, error, maximum. These three values are stored in private attributes 
	    and used to optimize the PS generator (basically use information from a successful run).
	.
	If results are written in the mentioned file this method return \c false otherwise \c true. This information
	is necessery to whether the integrator must be initialised or not.
  */	
  /*!
	\fn HD_PS_Base::AddChannel(std::string,double=1.)
	\brief Set pointers right

	This method calls HD_Channel_Selector::GetChannel and creates a pointer which after setting the
	integration channel weight according to the weight factor given is used to finally add this 
	(integration!) channel to the vector
	of channels in PHASIC::Multi_Channel class (PHASIC::Multi_Channel::Add).
  */
  /*!
	\fn void HD_PS_Base::Initialise() 
	\brief Initialisation of the PS integrator

	This method uses HD_PS_Base::CalculateNormalisedWidth and HD_PS_Base::WriteOut to 
	do that job.
  */	
  /*!
	\fn void HD_PS_Base::CalculateNormalisedWidth() 
	\brief Calculates optimised weight factors for Multi Channel integration

	In this method the weight factors of the integration channels are optimzed
	such that peak pattern in amplitudes are integrated with high efficiency.
  */	
  /*!
	\fn HD_PS_Base::WriteOut()
	\brief Writes integration channels into complete decay channel files.

	This method genererates the file corresponding to the considered decay channel.
	The information in this file, for instance, is picked out 
	using HD_PS_Base::Construct method in later runs.
  */	
  
  
} // end of namespace


#endif
