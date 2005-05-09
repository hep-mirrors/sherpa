#ifndef Three_Body_PSs_h
#define Three_Body_PSs_h

#include "Single_Channel.H"
#include "Vector.H"

namespace HADRONS {
  class Dalitz : public PHASIC::Single_Channel {
    ATOOLS::Vec4D m_decvec, m_pvec;
    double        m_pmass, m_pwidth, m_sexp, m_smin, m_smax;
    int           m_p1, m_p2, m_dir, m_mode;
    std::string   m_propmass;				
    int           m_chnumber;			
  public :
    Dalitz(const ATOOLS::Flavour * fl,const ATOOLS::Flavour prop,
		   const int p1, const int p2 );
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL,double * _ran=NULL);
    void GeneratePoint(ATOOLS::Vec4D * p,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL);
    void GenerateWeight(ATOOLS::Vec4D * p);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return name;}
  }; // end of class

  /*!
	\file Three_Body_PSs.H
	\brief Declares the class HADRONS::Dalitz

	This file can be found in directory <tt>PS_Library</tt>.
  */	
  /*!
	\class Dalitz
	\brief Tool to handle 3-particle PS

      This channel is meant to be the default 2-body integrator.
      There is no specific angular dependence, i.e. cos(theta) in the rest 
      frame is distributed uniformly.

	  This clas is a subclass of PHASIC::Single_Channel.

	  However, its member methods call PHASIC::Channel_Elements methods in order to deliver a point in
	  PS or the PS element (weight factor)
  */	  
  /*!
	\fn Dalitz::Dalitz(const ATOOLS::Flavour *,const ATOOLS::Flavour,const int,const int)
	\brief Constructor for Dalitz 1 to 3 decay integration channel

	The required information is in which two of the three final particles the
	propagator decays \f$(i,j)\f$ and of what kind the propagator is. 
  */	
  /*!
	\var Dalitz::m_decvec
	This is the 4-momentum of the decaying particle. Initially, its value is \f$(M,0,0,0)\f$.
  */  
  /*!
	\var Dalitz::m_pvec
	This is the 4-momentum of the propagator.
  */
  /*!
	\var Dalitz::m_pmass
	This is the mass of the propagator.
  */
  /*!
	\var Dalitz::m_pwidth
	This is the width of the propagator.
  */	
  /*!
	\var Dalitz::m_smin
	Minimum value \f$s_{min}=(m_i+m_j)^2\f$.
  */
  /*!
	\var Dalitz::m_smax
	Maximum value \f$s_{max}=(M-m_k)^2\f$ where \f$k \ne i, k \ne j\f$.
  */	
  /*!
	\var Dalitz::m_p1
	This corresponds to \f$i\f$.
  */	
  /*!
	\var Dalitz::m_p2
	This corresponds to \f$j\f$.
  */	
  /*!
	\var Dalitz::m_dir
	This corresponds to \f$k\f$.
  */	
  /*!
	\var Dalitz::m_mode
	Its value is 0 if the propagator mass is small enough, i.e.
	\f$m_p<\frac{1}{10}(m_i+m_j)\f$, which corresponds to a nearly
	massless propagator. Otherwise, it takes value 1.
  */	
	
	
  /*! 
	\fn Dalitz::ChID()
	Returns value in PHASIC::Single_Channel::name containing the keyword \c "Dalitz"
	followed by information about propagator.
  */	


} // end of namespace

#endif
