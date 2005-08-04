#ifndef HD_ME_Base_H
#define HD_ME_Base_H

#include "Flavour.H"
#include "Vector.H"
#include <string>

namespace HADRONS {

  struct GeneralParam {				// parameter
	double a, b;									// V-A structure
	double a2, b2;									// V-A structure
	double fpi, fK; 								// decay constants
	double Vud, Vus; 								// CKM ME; 
	double GF;										// Fermi constant
	double frho, grpp;								// rho-pipi coupling
	double alpha, beta, delta, gamma;				// weight factors for different resonances
	double Beta[2];									// weight factor for resonances
	double MR, GR, MRR, GRR, MO, GO, MRRR, GRRR;	// rho, omega resonances
	double MV[2], GV[2], Mv[2], Gv[2];				// vector resonances
	double MA, GA; 									// axial resonance a1(1260)
	double lambda;									// parameter in KS model
	double lambda0, exp_alpha;						// parameters in RChT parameterisation
	double gammaR, gammaRR, gammaO, gammaRRR;		// global factor for RChT off-shell width
  }; 

  struct GeneralVector {
	int i, j;
  };

  struct GeneralModel {
	int me;							// model for ME
	int run;						// running coupling
	int ff;							// form factor
	GeneralParam pm;			 	// coupling parameters
	GeneralVector vec[2];			// in which particles do vector resonances decay
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
