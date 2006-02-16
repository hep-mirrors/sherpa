#ifndef Four_Body_PSs_h
#define Four_Body_PSs_h

#include "Single_Channel.H"
#include "Vector.H"
#include "Single_Channel.H"
#include "Run_Parameter.H"
#include "Channel_Elements.H"
#include "Vegas.H"
#include "ResonanceFlavour.H"

namespace HADRONS {
  class TwoResonances : public PHASIC::Single_Channel {
	ATOOLS::Vec4D     m_P;
    int           	  m_i, m_j, m_k, m_dir;
    int           	  m_chnumber;			
	ResonanceFlavour  m_prop1, m_prop2;
	ATOOLS::Flavour	* p_fl;

	ATOOLS::Info_Key m_kI_123_4,m_kI_12_3,m_kI_1_2;
	PHASIC::Vegas* p_vegas;

  public :
	TwoResonances(
		const ATOOLS::Flavour * fl,
		ResonanceFlavour prop1, 
		const int _k,
		ResonanceFlavour prop2,
		const int _i,
		const int _j );
	~TwoResonances() {
	  if (p_fl!=NULL) delete [] p_fl; p_fl = NULL; 
	}
    void GeneratePoint(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL,double * _ran=NULL);
    void GeneratePoint(ATOOLS::Vec4D * p,double * _ran=NULL);
    void GenerateWeight(ATOOLS::Vec4D * p,ATOOLS::Cut_Data * =NULL);
    void GenerateWeight(ATOOLS::Vec4D * p);
    int  ChNumber()                 { return m_chnumber;       }
    void SetChNumber(int _chnumber) { m_chnumber = _chnumber;  }
    std::string ChID()              { return name;}
	void   AddPoint(double value )	{ 
	  Single_Channel::AddPoint( value );
	  p_vegas->AddPoint(value, rans); 
	}
    void   Optimize()                { p_vegas->Optimize(); } 
    void   EndOptimize()             { p_vegas->EndOptimize(); } 
    void   WriteOut(std::string pId) { p_vegas->WriteOut(pId); } 
    void   ReadIn(std::string pId)   { p_vegas->ReadIn(pId); } 
	void   CreatePointerOnVegas() {
//	  p_vegas = new PHASIC::Vegas( rannum, 100, name );
	}
	void   DeletePointerOnVegas() {
//	  delete p_vegas;
	}
  }; // end of class

  /*!
	\file Four_Body_PSs.H
	\brief Declares the class HADRONS::TwoResonances

	This file can be found in the directory <tt>PS_Library</tt>.
  */	
  /*!
	\class TwoResonances
	\brief Tool to handle 4-particle PS with 2 resonances

	This class is a subclass of PHASIC::Single_Channel.
     
    It is the default integrator for decays of the form with two resonances (internal propagators)
    \f[ M \to l + res_1; \: res_1 \to k+res_2; \: res_2 \to i+j. \f]
  */	  
  /*!
	\fn TwoResonances::TwoResonances ( const ATOOLS::Flavour * fl, ResonanceFlavour prop1, const int _k, ResonanceFlavour prop2, const int _i, const int _j )
	\brief Constructor for the TwoResonances 1 to 4 decay integration channel

	The required information is the structure of the decay, i.e.
	\f$M \to l + res_1 \to l + k + res_2 \to l + k + i + j\f$.
  */	
  /*!
	\var TwoResonances::m_P
	This is the 4-momentum of the decaying particle. Initially, its value is \f$(M,0,0,0)\f$.
  */  
  /*!
	\var TwoResonances::m_i
	This corresponds to \f$i\f$.
  */	
  /*!
	\var TwoResonances::m_j
	This corresponds to \f$j\f$.
  */	
  /*!
	\var TwoResonances::m_k
	This corresponds to \f$k\f$.
  */	
  /*!
	\var TwoResonances::m_dir
	This corresponds to \f$l\f$.
  */	
	
	
  /*! 
	\fn TwoResonances::ChID()
	Returns value in PHASIC::Single_Channel::name containing the keyword \c "TwoResonances"
	followed by information about the propagators.
  */	


} // end of namespace

#endif
