#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"
#include "Histogram.H"
#include "Message.H"
#include "Tools.H"

namespace HADRONS {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Lepton : public HD_ME_Base {
	int m_nutau, m_nulep, m_lep;
	double GF;					    // Fermi constant
	Complex cR1, cL1, cR2, cL2;		// couplings for Z function
	double a, b, a2, b2;
	public:
	  Tau_Lepton( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p ); 
	  double Using_Traces_Simple( const ATOOLS::Vec4D *_p );
	  double Using_Traces( const ATOOLS::Vec4D *_p );
	  double Using_Traces_2( const ATOOLS::Vec4D *_p );
	  void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Pseudo : public HD_ME_Base {
	int m_nutau, m_pion;
	double fxx;					// decay constants (f_pi, f_K)
	double Vxx;					// CKM ME (Vud, Vus)
	double GF;					// Fermi constant
	Complex cR, cL;				// coupling for X function
	public:
	  Tau_Pseudo( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  inline void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Two_Pion : public HD_ME_Base {
	struct Model {
	  int run;						// running width of rho
	  int ff;						// form factor
	}; 
	int m_nutau, m_pion_ch, m_pion0;
	struct Model m_md;
	double Vud;					// CKM ME
	double fxx;					// decay constant
	double GF;					// Fermi constant
	double frho, grpp;			// rho-pipi coupling
	double CG;					// Clebsch-Gordon-coefficient
	Complex cR, cL;				// coupling for X function
	double MR, GR, MO, GO, MR2, MO2;
	double beta, lambda, delta, gamma;
	double MRR, GRR, MRRR, GRRR, MRR2, MRRR2;
	double gammaR, gammaRR, gammaRRR;
	double m, m2;				// mass of pion or kaon
	public:
	  Tau_Two_Pion( int nout, ATOOLS::Flavour *fl );
	  double operator()(const ATOOLS::Vec4D *_p );
	  Complex A( double x, double y );
	  Complex FormFactor( double s ); 
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Pion_Kaon : public HD_ME_Base {
	 
	class FF_Base {
	  public :
	  FF_Base() { };
	  virtual Complex VectorFormFactor( double ) = 0;
	  virtual Complex ScalarFormFactor( double ) = 0;
	  virtual void SetModelParameters( GeneralModel ) = 0;
	  virtual void SetMasses2( double, double, double ) {};
	};
	  

	class RChT : public FF_Base {		// Guerrero, Tesis Doctoral
	   
	  double MK2, GK, fpi2;             // K* resonance par's, pi decay constant
	  double MK02, GK0;                 // K0* resonance par's
	  double renorm2;                   // renormalisation scale
	  double mPi2, mK2, mEta2, mPi, mK, mEta;  // mass of pi, K, eta
	  double Sigma_KP, Delta_KP;               // sum, difference of masses^2
	  
	  Complex JBar( double s, double MP2, double MQ2, double Sigma, double Delta );
	  Complex JBarBar( double s, double MP2, double MQ2, double Sigma, double Delta );
	  Complex Mr( double s, double MP2, double MQ2 );
	  Complex L( double s, double MP2, double MQ2 );
	  double MassWidthVector( double s );
	  double MassWidthScalar( double s );
	   
	  public:
	  RChT() {};
	  ~RChT() {};
	  Complex VectorFormFactor( double );
	  Complex ScalarFormFactor( double );
	  void SetModelParameters( GeneralModel _md );
	  void SetMasses2( double, double, double );
	}; 

	class KS : public FF_Base {		// Guerrero, Tesis Doctoral
	   
	  public:
	  KS() {};
	  ~KS() {};
	  void SetModelParameters( GeneralModel _md );
	  Complex VectorFormFactor( double );
	  Complex ScalarFormFactor( double );
	}; 

	int m_nutau, m_pion, m_kaon;
	double Vus2;				// CKM ME
	double fpi2;				// decay constant
	double GF2;					// Fermi constant
	Complex cR, cL;				// coupling for X function
	double m_ms[4];             // masses^2 
	double Delta_KP;            // difference of masses^2

	FF_Base * p_ff;
	 
	public:
	  Tau_Pion_Kaon( int nout, ATOOLS::Flavour *fl );
	  ~Tau_Pion_Kaon() { if (p_ff!=NULL) delete p_ff; }
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Three_Pseudo : public HD_ME_Base {
	 
	class FF_Base {
	  public :
	  FF_Base() { };
	  virtual void SetPath( std::string ) = 0;
	  virtual Complex FormFactor( int, double, double, double ) = 0;
	  virtual void SetModelParameters( GeneralModel ) = 0;
	  virtual void SetMode( int ) {};
	  virtual void SetInternalNumbers( int * ) {};
	  virtual void SetOutgoingMasses2( double, double, double ) {};
	};
	  

	class RChT : public FF_Base {		// Dumm, Pich, Portoles, hep-ph/0312183 
	  inline double Tau_Three_Pseudo::RChT::Sqrt_Lambda( double _a, double _b, double _c ) {
		return sqrt(Tools::Lambda(_a,_b,_c)); }
 
	  double  MassWidthVector( double s ); 
	  double  MassWidthAxial( double Q2 ); 
	   
	  void    CreatePhiHistogram();		// pre-calculation of phi(q2)
	  double  Phi( double Q2 ); 
	  double  IntegralPhi( double );
	   
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
	  RChT() { };
	  ~RChT() { if (p_phi!=NULL) delete p_phi; }
	  Complex FormFactor( int, double, double, double ); 
	  void SetModelParameters( GeneralModel _md );
	  void SetPath( std::string s ) { m_path = s; }
	}; 

	class KS : public FF_Base {		// Decker, Finkemeier, Mirkes, hep-ph/9310270
	  // order of Pseudos: see table 1 of the paper
	  inline Complex Tau_Three_Pseudo::KS::BW_V( int a, int b, double s ) {
		return Tools::BreitWigner( s, msV[a],  widthV[a], ms[a], ms[b], 1. ); }
	  inline Complex Tau_Three_Pseudo::KS::BW_v( int a, int b, double s ) {
		return Tools::BreitWigner( s, msv[a],  widthv[a], ms[a], ms[b], 1. ); }
	  inline double Sqrt_Lambda( double _a, double _b, double _c ) { 
		return sqrt(Tools::Lambda(_a,_b,_c)); }
		 
	  Complex BW_A( double s );
	  Complex Tvector1( int a, int b, double x );
	  Complex Tvector2( int a, int b, double x );
	  Complex TSvector( int a, int b, int c, double Q2, double s, double t );
	  Complex Tgen( int a, int b, int c, double s, double t);
	   
	  void    CreateGHistogram();		// pre-calculation of G(q2)
	  double  G( double s );
	  double  IntegralG( double );
	   
	  double MA2, MAA2;				// axial resonances
	  double GA, GAA;
	  double alpha, Beta[2];		// relative strength of resonances
	  double ms[3], msV[2], msv[2];	// masses of particles, vector resonances
	  double widthV[2], widthv[2];
	  int    running_width;	
	  ATOOLS::Histogram * p_G;    	// pointer on G histogram
	  std::string m_path;			// path of Decaydata files
	  int *  m_part;				// internal numbers for particles

	  int    m_mode;         		// internal mode representation
	  bool   m_twoident;			// true if two vector resonances are identical
	  bool   m_G123;				// false if there is only one amplitude
	  double m_X123, m_ms123;
	   
	  public:
	  KS() { } 
	  ~KS() { if (p_G!=NULL) delete p_G; }
	  Complex FormFactor( int, double, double, double );
	  void SetMode( int );
	  inline void SetInternalNumbers( int * part ) { m_part = part; }
	  inline void SetOutgoingMasses2( double ms1, double ms2, double ms3 ) {
		ms[0] = ms1;						// mass of first pseudo
		ms[1] = ms2;						// mass of second pseudo
		ms[2] = ms3;						// mass of third pseudo
	  }
	  void SetModelParameters( GeneralModel _md );
	  void SetPath( std::string s ) { m_path = s; }
	}; 

	int    m_nutau, m_pseudo_3, m_pseudo_1, m_pseudo_2;
	int    m_part[4];
	double Vud,Vus;				// CKM ME
	double fpi;					// decay constant
	double GF;					// Fermi constant
	Complex cR, cL;				// X function couplings
	double m_ms[5];				// masses^2
	int    m_mode;				// code for decay mode
	double m_A123;				// coupling constant
	void   SetA123();
	FF_Base * p_ff;
	public:
	   
	  Tau_Three_Pseudo( int nout, ATOOLS::Flavour *fl );
	  ~Tau_Three_Pseudo() {
		if (p_ff!=NULL) delete p_ff;
	  }
	 
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  Complex FormFactor( int j, double Q2, double s, double t ); 
	  void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Four_Pion_3 : public HD_ME_Base {

	class LorenzBase {
	  public :
	  Complex * X;
	   
		ATOOLS::Vec4D   * p;
		ATOOLS::Vec4D     r[5];
		double            s[5];
		double            q2;
	   
	  LorenzBase( ) {};
	  virtual Complex operator()( int ) {};
	  virtual void SetModelParameters( GeneralModel ) {};
	  void SetPrivates( Complex *, ATOOLS::Vec4D * );
	};

	class KS : public LorenzBase { 	// CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
	   

	  double  fpi2;					// pion decay constant
	  double  Go3p, Frho, grop;		// coupling constants
	  double  mpi2, mpi02;			// mass^2 of pion, pion0
	  double  MR, MRR, MRRR;		// mass of rho resonances
	  double  MR2, MRR2, MRRR2;		// mass^2 of rho resonances
	  double  GR, GRR, GRRR;		// width of rho resonances
	  double  MO, MO2, GO;		    // omega resonance
	  double  MS, MS2, GS;			// sigm0 resonance
	  double  MF, MF2, GF;			// f0 resonance
	  double  MA, MA2, GA;			// axial resonance
	  double  beta, gamma;			// relative strength of vector resonances
	  double  sigma;				// relative strength of vector resonance
	  double  R[5];					// important factor
	  Complex OmegaPi(),  Beta_opi[4];
	  Complex AonePi(),   Beta_api[4];
	  Complex SigmaRho(), Beta_srh[4];
	  Complex FzeroRho(), Beta_frh[4];
	   
	  double  Dots( int, int );
	  Complex Fk( double, Complex * );
	  Complex Trho( double );
	  Complex TTrho( double );
	   
	  public : 
	   
	  KS() : LorenzBase() {};
	  Complex operator()( int );
	  void SetModelParameters( GeneralModel _md );
	}; 

	static const int  m_ncontrib = 4;	// # contributions to total hadr current
	int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
	double       m_ms[6];
	int          m_inter[6];			// internal numbering
	 
	double       Vud2,Vus2;				// CKM ME
	double       GF2;					// Fermi constant
	Complex      cR, cL;				// couplings in X functions
	LorenzBase * p_lorenz;				// pointer to current contributions
	Complex		 m_Alpha[m_ncontrib];	// weight of each contribution
	Complex		 m_SumAlpha;			// sum of all weights
	 
	ATOOLS::Vec4D	m_p[6];				// momenta
	Complex 	    m_X[5];				// value of an X function
	 
	public:
	 
	Tau_Four_Pion_3( int nout, ATOOLS::Flavour *fl );
	~Tau_Four_Pion_3() {
	  if (p_lorenz!=NULL) delete p_lorenz;
	}
   
	double operator()(const ATOOLS::Vec4D *_p );
	double Using_Hels( const ATOOLS::Vec4D *_p );
	void SetModelParameters( GeneralModel _md );

  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Four_Pion_1 : public HD_ME_Base {

	class LorenzBase {
	  public :
	  Complex * X;
	   
		ATOOLS::Vec4D   * p;
		ATOOLS::Vec4D     r[5];
		double            s[5];
		double            q2;
	   
	  LorenzBase( ) {};
	  virtual Complex operator()() {};
	  virtual void SetModelParameters( GeneralModel ) {};
	  void SetPrivates( Complex *, ATOOLS::Vec4D * );
	};

	class KS : public LorenzBase { 	// CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
	   

	  double  fpi2;					// pion decay constant
	  double  mpi2, mpi02;			// mass^2 of pion, pion0
	  double  MR, MRR, MRRR;		// mass of rho resonances
	  double  MR2, MRR2, MRRR2;		// mass^2 of rho resonances
	  double  GR, GRR, GRRR;		// width of rho resonances
	  double  beta, gamma;			// relative strength of vector resonances
	   
	  double  Dots( int, int );
	  Complex Trho( double );
	   
	  public : 
	   
	  KS() : LorenzBase() {};
	  Complex operator()();
	  void SetModelParameters( GeneralModel _md );
	}; 

	static const int  m_ncontrib = 4;	// # contributions to total hadr current
	int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
	double       m_ms[6];
	int          m_inter[6];			// internal numbering
	 
	double       Vud2,Vus2;				// CKM ME
	double       GF2;					// Fermi constant
	Complex      cR, cL;				// couplings in X functions
	LorenzBase * p_lorenz;				// pointer to current contributions
	 
	ATOOLS::Vec4D	m_p[6];				// momenta
	Complex 	    m_X[5];				// value of an X function
	 
	public:
	 
	Tau_Four_Pion_1( int nout, ATOOLS::Flavour *fl );
	~Tau_Four_Pion_1() {
	  if (p_lorenz!=NULL) delete p_lorenz;
	}
   
	double operator()(const ATOOLS::Vec4D *_p );
	double Using_Hels( const ATOOLS::Vec4D *_p );
	void SetModelParameters( GeneralModel _md );

  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Eta_Two_Pion : public HD_ME_Base {
	 
	int    m_nutau, m_eta, m_pion_1, m_pion_2;
	double a, b;				// V-A structure
	double Vud;    				// CKM ME
	double fpi;					// decay constant
	double GF;					// Fermi constant
	Complex cR, cL;				// X function couplings
	double m_ms[5];				// masses^2
	 
	public:
	   
	  Tau_Eta_Two_Pion( int nout, ATOOLS::Flavour *fl );
	 
	  double operator()(const ATOOLS::Vec4D *_p );
	  double Using_Hels( const ATOOLS::Vec4D *_p );
	  void   SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Doxygen part                                                     %%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  /*!
	\class Tau_Three_Pseudo
	\brief Tools to calculate the ME for \f$\tau \to \nu_\tau \Pi \Pi \Pi\f$
   */	
  /*!
	\class Tau_Four_Pion
	\brief Tools to calculate the ME for \f$\tau \to \nu_\tau 4 \pi \f$
   */	
}
