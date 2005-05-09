#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"

namespace HADRONS {

  enum Models {
	MD_TRACES			= 1,
	MD_XYZ				= 2,
	MD_TRACES_SIMPLE	= 3
  }; 

  class Tau_Lepton : public HD_ME_Base {
	int m_nutau, m_nulep, m_lep;
	public:
	  Tau_Lepton( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Traces( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p ); 
  };

  class Tau_Pion : public HD_ME_Base {
	int m_nutau, m_pion;
	public:
	  Tau_Pion( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces( const ATOOLS::Vec4D *_p ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p );
  };

  class Tau_Two_Pion : public HD_ME_Base {
	int m_nutau, m_pion_ch, m_pion0;
	public:
	  Tau_Two_Pion( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  Complex A( double x, double y );
	  Complex FormFactor( double s ); 
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces( const ATOOLS::Vec4D *_p ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p );
  };

  /*!
	\file Tau_Decays_MEs.H
	\brief Declares the classes HADRONS::Tau_Lepton, HADRONS::Tau_Pion, and HADRONS::Tau_Two_Pion

	This file can be found in directory \c ME_Library
   */

  /*!
	\class Tau_Lepton
	\brief Tools to calculate the ME for leptonic decay channel
   */	
  /*!
	\class Tau_Pion
	\brief Tools to calculate the ME for \f$\tau \to \nu_\tau \pi/K\f$
   */
  /*!
	\class Tau_Two_Pion
	\brief Tools to calculate the ME for \f$\tau \to \nu_\tau \pi \pi^0\f$
   */	

}
