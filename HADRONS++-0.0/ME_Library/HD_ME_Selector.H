#ifndef HD_Selector_H
#define HD_Selector_H


#include "HD_ME_Base.H"
#include "Vector.H"
#include <string>

namespace HADRONS {

  class HD_ME_Selector {

    int  m_mode;
	void HD_ME_Selector::SelectTauDecay( int nout, ATOOLS::Flavour * flavs, HD_ME_Base *& hdme ); 
    void SelectLightPseudoScalarDecay( int nout, ATOOLS::Flavour *flavs, HD_ME_Base *& hdme );
    void SetVector_For_2PS(int,int,int,ATOOLS::Flavour *,HD_ME_Base *);
    void SetDalitz_For_3P(int,ATOOLS::Flavour *,HD_ME_Base *);

    bool IsPseudoScalar(ATOOLS::Flavour &);
	
  public:

    HD_ME_Selector()      : m_mode(0) {};
    HD_ME_Selector(int m) : m_mode(m) {};
    HD_ME_Base * GetME(int,int,ATOOLS::Flavour *,
		       std::string _met );
  }; // end of class

  inline bool HD_ME_Selector::IsPseudoScalar(ATOOLS::Flavour & flav) {
    int kfc = flav.Kfcode(); 
    if (kfc<1000 && kfc>110 && kfc%10==1) return true;
    return false;
  };  // end of inline (belongs to class)

  /*!
	\file HD_ME_Selector.H
	\brief Declares the class HADRONS::HD_ME_Selector

	This file can be found in directory \c ME_Library
  */
  /*!
	\class HD_ME_Selector
	\brief Tool to select the appropriate ME
  */	
  /*!
	\fn HD_ME_Selector::HD_ME_Selector()
	\brief Constructor with mode=0
  */
  /*!
	\fn HD_ME_Selector::HD_ME_Selector(int m)
	\brief Constructor with mode=m
  */
  /*!
	\var HD_ME_Selector::m_mode
	Mode number.
  */	
  /*!
	\fn HD_ME_Selector::GetME(int,int,ATOOLS::Flavour *,std::string _met=std::string(""))
	\brief Returns pointer on the corresponding ME.

	Firstly, the methods checks if the sum of outgoing masses is less than the 
	incoming masses. If not it returns NULL. Otherwise it returns a pointer on
	the corresponding ME class: 
	- Pseudoscalar -> 2 photons: P_2Gamma
	- Pseudoscalar -> 2 pseudoscalars + photon: P_2PGamma
	- Pseudoscalar -> 3 pseudoscalars: P_3P_Dalitz
	- Pseudoscalar -> photon + lepton + anti: P_GammaFF
	- Pseudoscalar -> pseudoscalar + 2 photons: P_P2Gamma
	.
	If none of the above list applies: the methods returns a pointer on Isotropic.

	This method uses HD_ME_Selector::SelectLightPseudoScalarDecay 
	method on order to succeed.
  */	
  /*!
	\fn HD_ME_Selector::SelectLightPseudoScalarDecay(int,ATOOLS::Flavour *,HD_ME_Base *&)
	Selects the corresponding ME class depending on given information (number of
	outgoing particles, flavours, etc). The pointer is mapped to the reference 
	argument
	<var>hdme</var> so that it is available in HD_ME_Selector::GetME 
	method which returns it
	to the Hadron_Decay_Channel constructor.
  */	




} // end of namespace

#endif
