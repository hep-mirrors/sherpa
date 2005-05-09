#ifndef Two_Body_MEs_H
#define Two_Body_MEs_H

#include "HD_ME_Base.H"

namespace HADRONS {
  class P_2Gamma : public HD_ME_Base {
  public:
    P_2Gamma(int,ATOOLS::Flavour *);
    double operator()(const ATOOLS::Vec4D *);
  }; // end of class

  /*!
	\file Two_Body_MEs.H
	\brief Declares the class HADRONS::P_2Gamma

	This code can be found in directory \c ME_Library.
  */

  /*!
	\class P_2Gamma 
	\brief Tool to calculate ME for \f$Pseudoscalar \to 2\gamma\f$ decay.

	This class is subclass of HADRONS::HD_ME_Base.
  */

  /*!
	\fn P_2Gamma::operator()(const ATOOLS::Vec4D *)
	\brief Returns amplitude squared.

	Up to now it does nothing else than to return \f$\frac{M}{4}\f$.
  */	

	
} // end of namespace


#endif
