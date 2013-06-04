#ifndef HADRONS_PS_Library_Three_Body_PSs_H
#define HADRONS_PS_Library_Three_Body_PSs_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Math/Vector.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"

namespace HADRONS {
  class Dalitz : public PHASIC::Single_Channel {

    private:
      ATOOLS::Vec4D m_decvec, m_pvec;
      double        m_pmass, m_pwidth, m_sexp, m_smin, m_smax;
      int           m_p1, m_p2, m_dir, m_mode;
      std::string   m_propmass;				
      int           m_chnumber;			
    public :
      Dalitz(
          const ATOOLS::Flavour * fl,
          SimpleResonanceFlavour res,
          const int p1, 
          const int p2 );
      void GeneratePoint(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL,
			 double * _ran=NULL);
      void GenerateWeight(ATOOLS::Vec4D * p,PHASIC::Cut_Data * =NULL);
      int  ChNumber()                 { return m_chnumber;       }
      void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
      std::string ChID()              { return name;}
  }; // end of class

  /*!
    \file Three_Body_PSs.H
    \brief Declares the class HADRONS::Dalitz

    This file can be found in the directory <tt>PS_Library</tt>.
  */	
  /*!
    \class Dalitz
    \brief Tool to handle 3-particle PS

    This class is a subclass of PHASIC::Single_Channel.
     
    It is the default integrator for decays of the form
    \f[ M \to k + res; \: res \to i+j. \f]
  */	  
  /*!
    \fn Dalitz::Dalitz( const ATOOLS::Flavour * fl, 
    SimpleResonanceFlavour res, const int p1, const int p2 )
    \brief Constructor for the Dalitz 1 to 3 decay integration channel
    
    The required information is in which two of the three final particles the
    propagator decays \f$(i,j)\f$ and of what kind the propagator is. 
  */	
  /*!
    \var Dalitz::m_decvec
    This is the 4-momentum of the decaying particle. Initially, its value is 
    \f$(M,0,0,0)\f$.
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
    Returns value in PHASIC::Single_Channel::name containing the keyword 
    \c "Dalitz"
    followed by information about propagator.
  */	


} // end of namespace

#endif
