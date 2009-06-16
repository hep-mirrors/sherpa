#ifndef AHADIC_Tools_Constituents_H
#define AHADIC_Tools_Constituents_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>


namespace AHADIC {


  class ConstituentCharacteristic {
  public:
    double m_mass;
    int    m_ispin;
    double m_weight, m_smearing;
    ConstituentCharacteristic(double _m=0., int _is=-1, 
			      double _wt=0., double _sm=0) : 
      m_mass(_m), m_ispin(_is), m_weight(_wt), m_smearing(_sm/_m) {}

    double Mass()      { return m_mass;     }
    double TotWeight() { return m_weight;   }
    double Smearing()  { return m_smearing; }
    int    ISpin()     { return m_ispin;    }
  };




  typedef std::map<ATOOLS::Flavour,ConstituentCharacteristic *> FlavCCMap;
  typedef FlavCCMap::iterator FlavCCMap_Iterator;

  class Constituents {
    double m_minmass;
  public:
    FlavCCMap CCMap;
    Constituents(bool=false);
    ~Constituents();

    double Mass(const ATOOLS::Flavour & flav);
    double MinMass();
    double TotWeight(const ATOOLS::Flavour & flav);
    double Smearing(ATOOLS::Flavour & flav);
    int    ISpin(ATOOLS::Flavour & flav);

    void   PrintConstituents();
  };


  /*!
    \file 
    \brief Contains the classes AHADIC::ConstituentCharacteristic and AHADIC::Constituents.
  */

  /*!
    \class ConstituentCharacteristic
    \brief Includes all information (constituent mass, spin d.o.f., flavour weight, 
    spin weight and angular smearing parameter) for a flavour that is needed in the cluster 
    hadronisation model. It is accessible through the class Constituents. 

    This class is initialised by the AHADIC::Hadronisation_Parameters and made available to 
    the rest of the code through a corresponding Get-method, 
    AHADIC::Hadronisation_Parameters::GetConstituents().

    The total weights are given in the following way:
    -  u, d quarks:   \f$wt = 1\f$ 
    -  s quarks: \f$wt = \chi_S\f$ 
    -  \f$ud_0\f$ diquarks: \f$wt = \chi_B\f$ 
    -  \f$sd_0\f$, \f$su_0\f$ diquarks: \f$wt = \chi_{SQ}\chi_B\f$ 
    -  \f$dd_1\f$, \f$uu_1\f$ diquarks: \f$wt = 3\chi_1\chi_B\f$ 
    -  \f$sd_1\f$, \f$su_1\f$ diquarks: \f$wt = 3\chi_1\chi_{SQ}\chi_B\f$ 
    -  \f$ss_1\f$ diquarks: \f$wt = 3\chi_1\chi_{SS}\chi_B\f$ 

    Here \f$\chi_{B, S}\f$ are the baryon and strangeness suppression parameters, respectively.
    With the suppression factors \f$\chi_{SQ}\f$ and \f$\chi_{SS}\f$ for dingle and double
    strange diquarks and \f$chi_1\f$ for spin-1 diquarks, the sum of weights for the different
    diquarks species \f$\Sigma\f$ is given by

    \f$\Sigma = 2+\chi_S+\chi_B\cdot[1+2\chi_{SQ}+3\chi_1(2+2\chi_{SQ}+\chi_{SS})]\f$.

    Since the methods are just simple access methods, no further documentation seems to be 
    neccessary.
  */
  /*!
    \class Constituents
    \brief A wrapper around a searchable map, where flavours are connected with their respective 
    ConstituentCharacteristic.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Constituents
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var FlavCCMap Constituents::CCMap
    \brief The map of flavours and ConstituentCharacteristic's, on which this class operates and
    which is its actual core.
  */
  /*!
    \var double m_minmass
    \brief The mass of the lightest constituent (no gluon), usually a d or u quark.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Constituents
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Constituents::Constituents(bool=false)
    \brief Constructor, the flag indicates whether diquarks are included. 

    The Constituents are constructed by the Hadronisation_Parameters, and actually use them in 
    order to fill the Constituents::CCMap. To this end, an explicit iteration is made over all 
    (light) flavours where masses and parts of the flavour weights are taken from the 
    Hadronisation_Parameters and the spin weights are fixed in the code (2 for quarks, 1 or 3 
    for scalar or vector diquarks, respectively).
  */
  /*!
    \fn Constituents::~Constituents()
    \brief Deletes the CCMap.
  */
  /*!
    \fn double Constituents::MinMass()
    \brief The constituent mass of the lightest constituent.
  */
  /*!
    \fn double Constituents::Mass(ATOOLS::Flavour & flav)
    \brief The constituent mass of a flavour.

    If the flavour is found in the CCMap, its constituent mass is returned; 
    if not, its original mass from the class ATOOLS::Flavour is returned.
  */
  /*!
    \fn double Constituents::FlWeight(ATOOLS::Flavour & flav)
    \brief The flavour weight for a given flavour.

    If the flavour is found in the CCMap, its flavour weight is returned; if not the 
    method yields 0.
  */
  /*!
    \fn double Constituents::SpWeight(ATOOLS::Flavour & flav)
    \brief The spin weight of a flavour.

    If the flavour is found in the CCMap, its spin weight, i.e. the number of spin degrees of 
    freedom, is returned; if not the method yields 0.
  */
  /*!
    \fn double Constituents::TotWeight(ATOOLS::Flavour & flav)
    \brief The total weight for a given flavour.

    If the flavour is found in the CCMap, its total weight (namely the total weight, given in
    AHADIC::ConstituentCharacteristic = flavour weight times spin weight) is returned; if not 
    the method yields 0.
  */
  /*!
    \fn double Constituents::Smearing(ATOOLS::Flavour & flav)
    \brief The smearing parameter for a given flavour.

    If the flavour is found in the CCMap, its angular smearing parameter is returned; if not 
    the method yields 0.  This is not in use add the moment.  
  */
  /*!
    \fn double Constituents::ISpin(ATOOLS::Flavour & flav)
    \brief The integer spin (twice the physical spin) for a given flavour.

    If the flavour is found in the CCMap, its integer spin is returned; if not the method 
    yields 0.
  */
  /*!
    \fn double Constituents::PrintConstituents()
    \brief Prints the ingedients of the CCMap.
  */
}

#endif
