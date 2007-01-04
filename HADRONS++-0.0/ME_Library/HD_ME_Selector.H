#ifndef HD_Selector_H
#define HD_Selector_H


#include "HD_ME_Base.H"
#include "Vector.H"
#include <string>

namespace HADRONS {

  class HD_ME_Selector {

    void SelectTauDecay(int nout,ATOOLS::Flavour * flavs, HD_ME_Base *& hdme ); 
    void SelectLightPseudoScalarDecay( int nout, ATOOLS::Flavour *flavs, HD_ME_Base *& hdme );
    void SetVector_For_2PS(int,int,int,ATOOLS::Flavour *,HD_ME_Base *);

    bool IsPseudoScalar(ATOOLS::Flavour &);
	
  public:

    HD_ME_Selector() {};
    HD_ME_Base * GetME(int nin,int nout,ATOOLS::Flavour * flavs );
  }; // end of class

  inline bool HD_ME_Selector::IsPseudoScalar(ATOOLS::Flavour & flav) {
    int kfc = flav.Kfcode(); 
    if (kfc<1000 && kfc>110 && kfc%10==1) return true;
    return false;
  }  // end of inline (belongs to class)

  /*!
	\file HD_ME_Selector.H
	\brief Declares the class HADRONS::HD_ME_Selector

	This file can be found in the directory \c ME_Library
  */
  /*!
	\class HD_ME_Selector
	\brief Tool to select the appropriate ME
  */	
  /*!
	\fn HD_ME_Selector::HD_ME_Selector()
	\brief Constructor 
  */
  /*!
	\fn HD_ME_Selector::GetME(int nin,int nout,ATOOLS::Flavour * flavs)
	\brief Returns pointer on the corresponding ME.

	Initially, the methods checks if the sum of outgoing masses is less than the 
	incoming masses. If not it returns NULL. Otherwise it returns a pointer on
	the corresponding ME class. 

    If there is no ME class implemented yet,
	the methods returns a pointer on Isotropic.

	This method uses 
    - HD_ME_Selector::SelectLightPseudoScalarDecay 
    - HD_ME_Selector::SelectTauDecay
    - HD_ME_Selector::SelectBMesonDecay
    - HD_ME_Selector::SelectKMesonDecay
    - etc.
    .
    to find the corresponding ME class for the decay channel.
  */	
  /*!
	\fn HD_ME_Selector::SelectLightPseudoScalarDecay(int,ATOOLS::Flavour *,HD_ME_Base *&)
	Selects the corresponding ME class depending on given information (number of
	outgoing particles, flavours, etc). The pointer is mapped to 
	<var>hdme</var> so that it is available in the HD_ME_Selector::GetME 
	method which returns it
	to the Hadron_Decay_Channel constructor.
  */	
  /*!
	\fn HD_ME_Selector::SelectTauDecay(int,ATOOLS::Flavour *,HD_ME_Base *&)
	Selects the corresponding ME class depending on given information (number of
	outgoing particles, flavours, etc). The pointer is mapped to 
	<var>hdme</var> so that it is available in the HD_ME_Selector::GetME 
	method which returns it
	to the Hadron_Decay_Channel constructor.
  */	




} // end of namespace

#endif
