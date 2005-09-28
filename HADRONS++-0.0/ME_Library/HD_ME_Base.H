#ifndef HD_ME_Base_H
#define HD_ME_Base_H

#include "Flavour.H"
#include "Vector.H"
#include <string>

namespace HADRONS {

  typedef std::map<std::string,double> String_Map;

  struct GeneralParam {				// parameter
	double a, b;									// V-A structure
	double a2, b2;									// V-A structure
	double fpi, fK; 								// decay constants
	double Vud, Vus; 								// CKM ME; 
	double Go3p;      								// coupling constants for tau->4pion
	double GF;										// Fermi constant
	double frho, grpp;								// rho-pipi coupling
	double grop;									// rho-omega-pi coupling
	double alpha, beta, delta, gamma;				// weight factors for different resonances
	double sigma;									
	double Beta[2];									// weight factor for resonances
	Complex Alpha[4];					
	double MR, GR, MRR, GRR, MO, GO, MRRR, GRRR;	// rho, omega resonances 
	double MK, GK;                                  // K* resonance
	double MK0, GK0;                                // K0* resonance
	double MV[2], GV[2], Mv[2], Gv[2];				// vector resonances
	double MA, GA, MAA, GAA;						// axial resonance 
	double Mdummy2, Mdummy1;						// mass of (e.g.) sigma, f0
	double Gdummy2, Gdummy1;						// width of (e.g.) sigma, f0 
	double dummy1;									// dummy for whatsoever
	 
	double lambda;									// parameter in KS model
	double lambda0, exp_alpha;						// parameters in RChT parameterisation
	double gammaR, gammaRR, gammaO, gammaRRR;		// global factor for RChT off-shell width
	Complex Beta_opi[4], Beta_api[4], Beta_frh[4], Beta_srh[4]; // relative strength needed in tau->4pion
  }; 

  struct GeneralVector {			// in which two particles does vector resonance decay
	int i, j;
  };

  struct GeneralModel: public String_Map {
	int run;						// running coupling
	int ff;							// form factor
	GeneralParam pm;			 	// coupling parameters
	GeneralVector vec[2];			// in which particles do vector resonances decay
  public:
    inline double operator()(const std::string &tag,const double &def) const 
    {
      String_Map::const_iterator fit(find(tag));
      return fit!=end()?fit->second:def;
    }
  }; 

  class HD_ME_Base {
  protected:
    std::string       m_metype;
    int               m_nout;
    ATOOLS::Flavour * p_flavs;
    double          * p_masses,* p_masses2;
	std::string       m_path;
  public:
    HD_ME_Base(int,ATOOLS::Flavour *);
    virtual ~HD_ME_Base();
    void SetPath( std::string s ) { m_path=s; }
	std::string METype() { return m_metype; }

    virtual double operator()(const ATOOLS::Vec4D *) { return 1.; }
    virtual void AddVector(double,double) {}
	virtual void SetModelParameters( struct GeneralModel _md ) {};
  }; // end of class HD_ME_Base

  class Isotropic : public HD_ME_Base {
  public:
    Isotropic(int,ATOOLS::Flavour *,
	      std::string _met=std::string(""));
  }; // end of class Isotropic
  
  /*!
	\file HD_ME_Base.H
	\brief Declares the classes HADRONS::HD_ME_Base and HADRONS::Isotropic

	This file can be found in directory \c ME_Library.
  */

  /*!
	\class HD_ME_Base
	\brief Tool to set up the connection to the corresponding ME method
  */	
  /*!
	\fn HD_ME_Base::HD_ME_Base
	\brief Constructor and initialisation of private attributes
  */	
  /*!
	\var HD_ME_Base::p_masses
	Pointer of masses.
  */	
  /*!
	\var HD_ME_Base::p_masses2
	Pointer of masses squared.
  */	
  /*!
	\var HD_ME_Base::p_flavs
	Pointer on flavours.
  */
  /*!
	\var HD_ME_Base::m_nout
	Number of outgoing particles.
  */
  /*!
	\var HD_ME_Base::m_metype
	Type of ME.
  */
  /*!
	\var HD_ME_Base::m_md
	Contains information which model has to be used for calculating a matrix element (only valid for tau decays).
	- Matrix element:
		- traces: polarization of tau by \f$S_\tau=(0,0,0,1)\f$
		- simple traces: no polarization
		- XYZ functions: no polarization
		.
	- Running width: the width of intermediate resonance particles is either fixed or depends on the momentum transfer
	- Parameter: V-A structur for electroweak currents (tau and other current in case of leptonic decay seperately).
	.
   */

  /*!
	\class Isotropic
	\brief Not worth mentioning

	This class is a subclass of HADRONS::HD_ME_Base.
  */	

} // end of namespace


#endif
