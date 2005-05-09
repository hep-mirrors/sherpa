using namespace ATOOLS;
using namespace std;

#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

namespace HADRONS
{

  class XYZFunc
  {
	  Vec4D  	* p_mom;
	  Flavour   * p_flav;
	  int N;
	  vector<Complex> eta, mu;

	  void CalcEtaMu() 
	  {
		Vec4D pi;
		mu.clear(); eta.clear();
		for( int i=0; i<N; i++ )
		{
		  pi = p_mom[i];
		  eta.push_back( csqrt( 2.*(pi[0]-SQRT_05*(pi[2]+pi[3])) ) );
		  Complex help( p_flav[i].Mass(), 0. );
		  mu.push_back( help/eta[i] );
		  if( p_flav[i].IsAnti() ) mu[i] *= -1.;
		}
	  }
	  
	public:

	  XYZFunc( int nout, const Vec4D *p, const Flavour *fl )
	  {
		N = nout+1;
		p_mom = new Vec4D[N];
		p_flav = new Flavour[N];
		for( int i=0; i<N; i++ ) {
		  p_mom[i] = p[i];
		  p_flav[i] = fl[i];
		}
		CalcEtaMu();
	  }

	  void Delete()
	  {
		delete [] p_mom;
		delete [] p_flav;
	  }

	  Complex Z( const int t1, const int t2, const int t3, const int t4,
		  const int hel_comb,
		  const Complex cR1, const Complex cL1,
		  const Complex cR2, const Complex cL2 );
	  Complex Y( const int t1, const int t2, const int hel_comb, const Complex cL, const Complex cR ); 
	  Complex X( const int t1, const int t2, const int t3, 
		  const int hel_comb, const Complex cL, const Complex cR ); 
	  Complex S( const int s, const int i, const int j );
	  Complex Q(short hel); 

  };


  /*!
	\file XYZFuncs.H
	\brief Declares the class HADRONS::XYZFunc

	This file can be found in directory \c ME_Library.
   */  

  /*!
	\class XYZFunc
	\brief Tools to calculate X, Y, and Z functions

	This class contains everything that is necessery to calculate the value
	of an X, Y, or Z function, which can be used to calculate a decay matrix
	element.
	<b>Note!</b> Helicity combinations are code in a binary number such that 
	- ++++ = 0 
	- +++- = 1 
	- ++-+ = 2 
	- ...
	- ---- = 15 
	.
   */  

  /*!
	\fn XYZFunc::XYZFunc( int nout, const Vec4D *p, const Flavour *fl )
	\brief Constructor to set up the calculation

	This method saves the momenta and flavours so that they are available throughout the class. Furthermore, it
	calls XYZFunc::CalcEtaMu which is tasked with the calculation of \f$\eta_i\f$ and \f$\mu_i\f$ for each
	particle.
   */  
  /*!
	\fn XYZFunc::CalcEtaMu()
	\brief Calculates \f$\eta_i\f$ and \f$\mu_i\f$ for each particle

	Eta's and Mu's are calcuted and stored in private attributes. For Eta's \f$k_0\f$ is set to
	\f$k_0=(1,0,\frac{1}{\sqrt{2}},\frac{1}{\sqrt{2}})\f$.
   */

  /*!
	\var XYZFunc::p_mom

	Pointer on involved 4-momenta.
   */

  /*!
	\var XYZFunc::p_flav

	Pointer on involved flavours.
   */

  /*!
	\var XYZFunc::N

	Number of involved particles (in- and out-particles)
   */	

	
  
} 
