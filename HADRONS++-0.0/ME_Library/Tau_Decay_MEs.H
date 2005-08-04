#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"
#include "Histogram.H"
#include "Message.H"

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

  class Tau_Pseudo : public HD_ME_Base {
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
	  Tau_Pseudo( int nout, ATOOLS::Flavour *fl );
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
	double MR, MP, GR, MO, GO;
	double MR2, MP2, MO2, MRR2, MRRR2;
	double alpha, beta, lambda, delta, gamma;
	double MRR, GRR, MRRR, GRRR;
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

		MP = ATOOLS::Flavour(ATOOLS::kf::pi_plus).Mass();
		MR = _md.pm.MR;
		GR = _md.pm.GR;
		MRR = _md.pm.MRR;
		GRR = _md.pm.GRR;
		MRRR = _md.pm.MRRR;
		GRRR = _md.pm.GRRR;
		MO = _md.pm.MO;
		GO = _md.pm.GO;

		MR2 = MR*MR;
		MRR2 = MRR*MRR;
		MRRR2 = MRRR*MRRR;
		MP2 = MP*MP;
		MO2 = MO*MO;
		 
		alpha = _md.pm.alpha;
		beta  = _md.pm.beta;
		delta = _md.pm.delta;
		lambda = _md.pm.lambda;
		gamma = _md.pm.gamma;
	  }
  };

  class DPP;
  class DFM; 
  class Tau_Three_Pseudo : public HD_ME_Base {
	 
	struct Model {
	  int me;						// model for ME
	  int run;						// running width of rho
	  int ff;
	}; 

	class DPP {						// Dumm, Pich, Portoles
	  double  MassWidthVector( double s ); 
	  double  MassWidthAxial( double Q2 ); 
	  Complex BreitWignerRho( double s ); 
	   
	  double  Sqrt_Lambda(double,double,double);
	  void    CreatePhiHistogram();		// pre-calculation of phi(q2)
	  double  Phi( double Q2 ); 
	  double  IntegralPhi( double );
	  double  FitOrder7oddLowPhi(double);
	  double  FitOrder4HighPhi(double);
	  double  FitOrder4InterPhi( double );
	  double  OptimisedPhi(double);
	   
	  double l0, exp_alpha;			// fitting parameter
	  double MV, MV2, MA, MA2;		// resonances (vector, axial)
	  double GA_at_MA2;				// on-shell axial width
	  double Phi_at_MA2;			// Phi(MA2)
	  double m, m2, mK2;			// masses of pion and kaon
	  double fpi;					// decay constant
	  double gammaR;				// global constant for GV
	  ATOOLS::Histogram * p_phi;    // pointer on phi histogram
	  std::string m_path;			// path of Decaydata files

	  public:
	  DPP() { };
	  ~DPP() { if (p_phi!=NULL) delete p_phi; }
	  Complex FormFactor( int, double, double, double ); 
	  void SetModelParameters( struct GeneralModel _md ) { 
		fpi  = _md.pm.fpi;
		l0   = _md.pm.lambda0;		   		// fit parameter lambda0
		MA   = _md.pm.MA;					// mass of axial resonance
		MV   = _md.pm.MV[0];				// mass of vector resonance
		MV2  = ATOOLS::sqr(MV);
		MA2  = ATOOLS::sqr(MA);
		gammaR = _md.pm.gammaR;
		m    = ATOOLS::Flavour( ATOOLS::kf::pi_plus ).Mass();
		m2   = ATOOLS::sqr(m);
		mK2  = ATOOLS::sqr( ATOOLS::Flavour( ATOOLS::kf::K_plus ).Mass() );
		exp_alpha = _md.pm.exp_alpha;		// exponent in off-shell GA
		GA_at_MA2 = _md.pm.GA;				// on-shell axial width
		CreatePhiHistogram();
		Phi_at_MA2 = Phi( MA2 );
	  }
	  void SetPath( std::string s ) { m_path = s; }
	}; 

	class DFM {					// Decker, Finkemeier, Mirkes
	  // order of Pseudos: see table 1 of the paper
	  Complex BW( double s, double M2, double MG );
	  Complex BW_V( int a, int b, double s );
	  Complex BW_v( int a, int b, double s );
	  Complex BW_A( double s );
	  Complex Tvector1( int a, int b, double x );
	  Complex Tvector2( int a, int b, double x );
	  Complex TSvector( int a, int b, int c, double Q2, double s, double t );
	  Complex Tgen( int a, int b, int c, double s, double t);
	  Complex Trho(double);
	   
	  void    CreateGHistogram();		// pre-calculation of G(q2)
	  double  G( double s );
	  double  IntegralG( double );
	  double  Sqrt_Lambda(double,double,double);
	   
	  double MA2;
	  double GA;
	  double Beta[2];
	  double ms[3], msV[2], msv[2];
	  double widthV[2], widthv[2];
	  int    running_width;	
	  ATOOLS::Histogram * p_G;    	// pointer on G histogram
	  std::string m_path;			// path of Decaydata files
	  int *  m_part;				// internal numbers for particles

	  int    m_mode;				// internal mode representation
	  bool   m_twoident;			// true if two vector resonances are identical
	   
	  public:
	  DFM() { } 
	  ~DFM() { if (p_G!=NULL) delete p_G; }
	  Complex FormFactor( int, double, double, double );
	  inline void SetOutgoingMasses2( double ms1, double ms2, double ms3 ) {
		ms[0] = ms1;						// mass of first pseudo
		ms[1] = ms2;						// mass of decons pseudo
		ms[2] = ms3;						// mass of third pseudo
	  }
	  inline void SetMode( int m ) { 
		m_mode = m;
		m_twoident = false;
		if (m_mode==1200 ||
		    m_mode==210 ||
			m_mode==30 ||
			m_mode==3000 ||
			m_mode==12) m_twoident = true;
	  }
	  inline void SetInternalNumbers( int * part ) { m_part = part; }
	  void SetModelParameters( GeneralModel _md ) {
		 
		int a[2], b[2], c[2];
		if (!m_twoident) {
		  for( int i=0; i<2; i++ ) {
			a[i] = m_part[ _md.vec[i].i-1 ];
			b[i] = m_part[ _md.vec[i].j-1 ];
			if ( a[i]!=3 && b[i]!=3 ) {
			  ATOOLS::msg.Error()<<ATOOLS::om::red
				<<"ERROR in Tau_Three_Pseudo::DFM::SetModelParameters\n"
				<<"     Resonances aren't set correctly.\n"
				<<"     Make sure you have the right settings under \"Resonances\" vector"<<i+1<<" -> _ _.\n"
				<<"     mode number : "<<m_mode<<"\n"
				<<"     Don't know what to, will abort."
				<<ATOOLS::om::reset<<std::endl;
			  abort();
			}
			else c[i] = (a[i]==3)? b[i] : a[i];
		  }
		}
		else {
		  c[0] = c[1] = 1;
		}
		MA2    = ATOOLS::sqr(_md.pm.MA);			// mass^2 of axial resonance
		msV[0] = ATOOLS::sqr(_md.pm.MV[c[0]-1]);	// mass^2 of vector resonance 13
		msV[1] = ATOOLS::sqr(_md.pm.MV[c[1]-1]);	// mass^2 of vector resonance 23
		msv[0] = ATOOLS::sqr(_md.pm.Mv[c[0]-1]);	// mass^2 of vector resonance' 13
		msv[1] = ATOOLS::sqr(_md.pm.Mv[c[1]-1]);	// mass^2 of vector resonance' 23
		Beta[0] = _md.pm.Beta[c[0]-1];				// weight factor for vector resonance' 13
		Beta[1] = _md.pm.Beta[c[1]-1];				// weight factor for vector resonance' 23
		 
		GA      = _md.pm.GA;						// on-shell axial width
		widthV[0] = _md.pm.GV[c[0]-1];				// on-shell vector 13 width
		widthV[1] = _md.pm.GV[c[1]-1];				// on-shell vector 23 width
		widthv[0] = _md.pm.Gv[c[0]-1];				// on-shell vector' 13 width
		widthv[1] = _md.pm.Gv[c[1]-1];				// on-shell vector' 23 width

		running_width = _md.run;

		CreateGHistogram();
	  } 
	  void SetPath( std::string s ) { m_path = s; }
	}; 

	int    m_nutau, m_pseudo_3, m_pseudo_1, m_pseudo_2;
	int    m_part[4];
	DPP    dpp;
	DFM    dfm;
	struct Model m_md;			
	double a, b;				// V-A structure
	double Vud,Vus;					// CKM ME
	double fpi;					// decay constant
	double GF;					// Fermi constant
	Complex cR, cL;				// X function couplings
	double m_ms[5];				// masses^2
	int    m_mode;				// code for decay mode
	double m_A123;				// coupling constant
	void   SetA123();
	public:
	  Tau_Three_Pseudo( int nout, ATOOLS::Flavour *fl );
	 
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  Complex FormFactor( int j, double Q2, double s, double t ); 
	  inline void SetModelParameters( struct GeneralModel _md ) { 
		m_md.me  = _md.me;
		m_md.run = _md.run;
		m_md.ff  = _md.ff;  
		a        = _md.pm.a;
		b        = _md.pm.b;
		Vud      = _md.pm.Vud;
		Vus      = _md.pm.Vus;
		fpi      = _md.pm.fpi;
		GF       = _md.pm.GF;
		cR = (0.,a - b);
		cL = (0.,a + b);
		SetA123();
		switch( m_md.ff ) {
		  case 2 : dpp.SetPath( m_path );
				   dpp.SetModelParameters( _md );
				   break;
		  case 1 : dfm.SetPath( m_path );
				   dfm.SetMode(m_mode);
		           dfm.SetOutgoingMasses2( m_ms[m_pseudo_1], 
		                                   m_ms[m_pseudo_2], 
										   m_ms[m_pseudo_3]  );
				   dfm.SetInternalNumbers( m_part );
				   dfm.SetModelParameters( _md );
				   break;
		}
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
