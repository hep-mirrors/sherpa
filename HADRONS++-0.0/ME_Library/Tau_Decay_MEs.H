#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"

namespace HADRONS {

  enum Models {
	MD_TRACES			= 1,
	MD_XYZ				= 2,
	MD_TRACES_SIMPLE	= 3
  }; 

  class Tau_Lepton : public HD_ME_Base {
	struct Model {
	  int me;
	}; 
	int m_nutau, m_nulep, m_lep;
	struct Model m_md;
	double a, b;					// V-A structure
	double a2, b2;				    // V-A structure
	double GF;					    // Fermi constant
	Complex cR1, cL1, cR2, cL2;		// couplings for Z function
	public:
	  Tau_Lepton( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Traces( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p ); 
	  inline void SetModelParameters( struct GeneralModel _md ) { 
		m_md.me  = _md.me;
		a    = _md.pm.a;
		b    = _md.pm.b;
		a2   = _md.pm.a2;
		b2   = _md.pm.b2;
		GF   = _md.pm.GF;
		cR1 = (0.,a-b);
		cL1 = (0.,a+b);
		cR2 = (0.,a2-b2);
		cL2 = (0.,a2+b2);
	  }
  };

  class Tau_Pion : public HD_ME_Base {
	struct Model {
	  int me;						// model for ME
	}; 
	int m_nutau, m_pion;
	struct Model m_md;
	double a, b;				// V-A structure
	double fxx;					// decay constants (f_pi, f_K)
	double Vxx;					// CKM ME (Vud, Vus)
	double GF;					// Fermi constant
	Complex cR, cL;				// coupling for X function
	public:
	  Tau_Pion( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces( const ATOOLS::Vec4D *_p ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  inline void SetModelParameters( struct GeneralModel _md ) { 
		m_md.me  = _md.me;
		a    = _md.pm.a;
		b    = _md.pm.b;
		Vxx  = ( p_flavs[m_pion].Kfcode() == ATOOLS::kf::pi_plus )? _md.pm.Vud : _md.pm.Vus;
		fxx  = ( p_flavs[m_pion].Kfcode() == ATOOLS::kf::pi_plus )? _md.pm.fpi : _md.pm.fK;
		GF   = _md.pm.GF;
		cR = (0.,a-b);
		cL = (0.,a+b);
	  }
  };

  class Tau_Two_Pion : public HD_ME_Base {
	struct Model {
	  int me;						// model for ME
	  int run;						// running width of rho
	  int ff;						// form factor
	}; 
	int m_nutau, m_pion_ch, m_pion0;
	struct Model m_md;
	double a, b;					// V-A structure
	double Vud;					// CKM ME
	double fxx;					// decay constant
	double GF;					// Fermi constant
	double frho, grpp;			// rho-pipi coupling
	double CG;					// Clebsch-Gordon-coefficient
	Complex cR, cL;				// coupling for X function
	public:
	  Tau_Two_Pion( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  Complex A( double x, double y );
	  Complex FormFactor( double s ); 
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces( const ATOOLS::Vec4D *_p ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  inline void SetModelParameters( struct GeneralModel _md ) { 
		m_md.me  = _md.me;
		m_md.run = _md.run;
		m_md.ff  = _md.ff; 
		a    	 = _md.pm.a;
		b    	 = _md.pm.b;
		Vud  	 = _md.pm.Vud; 
		fxx  	 = ( p_flavs[m_pion_ch].Kfcode() == ATOOLS::kf::pi_plus )? _md.pm.fpi : _md.pm.fK;
		GF   	 = _md.pm.GF;
		frho 	 = _md.pm.frho;
		grpp 	 = _md.pm.grpp;
		CG   	 = ( p_flavs[m_pion_ch].Kfcode() == ATOOLS::kf::pi_plus )? 1. : 0.5;
		cR       = (0.,a - b);
		cL       = (0.,a + b);
	  }
  };

  class Tau_Three_Pion : public HD_ME_Base {
	struct Model {
	  int me;						// model for ME
	  int run;						// running width of rho
	}; 
	int m_nutau, m_pion_ch, m_pion_1, m_pion_2;
	double  MassWidthVector( double s ); 
	double  MassWidthAxial( double Q2 ); 
	Complex BreitWignerRho( double s ); 
	double  Sqrt_Lambda( double a, double b, double c ); 
	double  Phi( double Q2 ); 
	double  IntegralPhi( double );
	double  FitOrder19oddPhi( double ); 
	double  FitOrder5Phi( double );
	Complex FormFactor( int j, double Q2, double s, double t, double u ); 
	double l0, exp_alpha;		// fitting parameter
	double MV, MV2, MA, MA2;	// resonances (vector, axial)
	double GA_at_MA2;			// on-shell axial width
	double Phi_at_MA2;			// Phi(MA2)
	double m, m2, mK2;			// masses of pion and kaon
	struct Model m_md;			
	double a, b;				// V-A structure
	double Vud;					// CKM ME
	double fpi;					// decay constant
	double GF;					// Fermi constant
	Complex cR, cL;				// X function couplings
	public:
	  Tau_Three_Pion( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  inline void SetModelParameters( struct GeneralModel _md ) { 
		m_md.me  = _md.me;
		m_md.run = _md.run;
		a        = _md.pm.a;
		b        = _md.pm.b;
		Vud      = _md.pm.Vud;
		fpi      = _md.pm.fpi;
		GF       = _md.pm.GF;
		MA   = 1.204;  					   // mass of axial resonance
		l0   = 11.9;					   // fit parameter lambda0
		MV   = 0.775;					   // mass of vector resonance
		MV2  = ATOOLS::sqr(MV);
		MA2  = ATOOLS::sqr(MA);
		m    = ATOOLS::Flavour( ATOOLS::kf::pi_plus ).Mass();
		m2   = ATOOLS::sqr(m);
		mK2  = ATOOLS::sqr( ATOOLS::Flavour( ATOOLS::kf::K_plus ).Mass() );
		exp_alpha = 2.45;
		GA_at_MA2 = 0.482799;
		Phi_at_MA2 = -10.3334;
		cR = (0.,a - b);
		cL = (0.,a + b);
	  }
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
