#ifndef Two_Body_PSs_h
#define Two_Body_PSs_h

#include "Single_Channel.H"
#include "Vector.H"

namespace HADRONS {
  class Iso2Channel : public PHASIC::Single_Channel {
    ATOOLS::Vec4D m_decvec;
    int           m_chnumber;
  public :
    Iso2Channel(const ATOOLS::Flavour *);
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL,double * _ran=NULL);
    void GeneratePoint(ATOOLS::Vec4D * p,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL);
    void GenerateWeight(ATOOLS::Vec4D * p);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return std::string("Iso2");}
  }; // end of class


  /*!
	\file Two_Body_PSs.H
	\brief Declares the class HADRONS::Iso2Channel
	
	This file can be found in directory \c PS_Library.
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
