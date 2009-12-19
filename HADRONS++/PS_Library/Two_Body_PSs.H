#ifndef HADRONS_PS_Library_Two_Body_PSs_H
#define HADRONS_PS_Library_Two_Body_PSs_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Math/Vector.H"

namespace HADRONS {
  class Iso2Channel : public PHASIC::Single_Channel {
    ATOOLS::Vec4D m_decvec;
    int           m_chnumber;
  public :
    Iso2Channel(const ATOOLS::Flavour *);
    void GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return std::string("Iso2");}
  }; // end of class

  class Iso1Channel : public PHASIC::Single_Channel {
    ATOOLS::Vec4D m_decvec;
    int           m_chnumber;
  public :
    Iso1Channel(const ATOOLS::Flavour *);
    void GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return std::string("Iso1");}
  }; // end of class


  /*!
	\file Two_Body_PSs.H
	\brief Declares the class HADRONS::Iso2Channel
	
	This file can be found in the directory \c PS_Library.
  */

  /*!
	\class Iso2Channel
	\brief Tool to handle isotropic 2-particle decay PS


      This channel is meant to be the default 2-body integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly. 
	  
	  This class is a subclass of PHASIC::Single_Channel

	  However, its member methods call PHASIC::Channel_Elements methods in order to deliver a point in
	  PS or the PS element (weight factor)
  */	
  /*!
	\fn Iso2Channel::Iso2Channel( const ATOOLS::Flavour * fl )
	\brief Constructor for isotropic 1 to 2 decay integration channel
  */
  /*!
	\var Iso2Channel::m_decvec
	This is the 4-momentum of the decaying particle. Initially, its value is \f$(M,0,0,0)\f$.
  */  
  /*! 
	\fn Iso2Channel::ChID()
	Returns value <tt>"Iso2"</tt>.
  */	
	
  

} // end of namespace

#endif
