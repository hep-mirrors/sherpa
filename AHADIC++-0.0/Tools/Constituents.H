#ifndef Constituents_H
#define Constituents_H

#include "Flavour.H"
#include <map>
#include <string>


namespace AHADIC {


  class ConstituentCharacteristic {
  public:
    double m_mass;
    int    m_ispin;
    double m_flweight, m_spweight, m_smearing;
    ConstituentCharacteristic(double _m=0., int _is=-1, 
			      double _fl=0., double _sp=0., double _sm=0) : 
      m_mass(_m), m_flweight(_fl), m_spweight(_sp), m_ispin(_is), 
      m_smearing(_sm/_m) {}

    double Mass()      { return m_mass;     }
    double FlWeight()  { return m_flweight; }
    double SpWeight()  { return m_spweight; }
    double TotWeight() { return m_spweight*m_flweight; }
    double Smearing()  { return m_smearing; }
    int    ISpin()     { return m_ispin;    }

    void   SetMass(double _mass)         { m_mass     = _mass;     }
    void   SetFlWeight(double _flweight) { m_flweight = _flweight; }
    void   SetSpWeight(double _spweight) { m_spweight = _spweight; }
    void   SetSmearing(double _smearing) { m_smearing = _smearing; }
    void   SetISpin(int _ispin)          { m_ispin    = _ispin;    }
  };




  typedef std::map<ATOOLS::Flavour,ConstituentCharacteristic *> FlavCCMap;
  typedef FlavCCMap::iterator FlavCCMap_Iterator;

  class Constituents {
  public:
    FlavCCMap CCMap;
    Constituents(bool=false);
    ~Constituents();

    double Mass(ATOOLS::Flavour & flav);
    double FlWeight(ATOOLS::Flavour & flav);
    double SpWeight(ATOOLS::Flavour & flav);
    double TotWeight(ATOOLS::Flavour & flav);
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
    \brief Includes all information (constituent mass, spin d.o.f., flavour weight, spin weight and 
    angular smearing parameter) for a flavour that is needed in the cluster hadronisation model. It is
    accessible through the class Constituents.

    Since the methods are just simple access methods, no further documentation seems to be neccessary.
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
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Constituents
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Constituents::Constituents(bool=false)
    \brief Constructor, the flag indicates whether diquarks are included. 

    The Constituents are constructed by the Hadronisation_Parameters, and actually use them in order to
    fill the Constituents::CCMap. To this end, an explicit iteration is made over all (light) flavours
    where masses and parts of the flavour weights are taken from the Hadronisation_Parameters and the spin
    wights are fixed in the code (2 for quakrs, 1 or 3 for scalar or vector diquarks, respectively).
  */
  /*!
    \fn Constituents::~Constituents()
    \brief Deletes the CCMap.
  */
  /*!
    \fn double Constituents::Mass(ATOOLS::Flavour & flav)
    \brief If the flavour is found in the CCMap, its constituent mass is returned; if not, its original mass
    from the class ATOOLS::Flavour is returned.
  */
  /*!
    \fn double Constituents::FlWeight(ATOOLS::Flavour & flav)
    \brief If the flavour is found in the CCMap, its flavour weight is returned; if not the method yields 0.
  */
  /*!
    \fn double Constituents::SpWeight(ATOOLS::Flavour & flav)
    \brief If the flavour is found in the CCMap, its spin weight is returned; if not the method yields 0.
  */
  /*!
    \fn double Constituents::TotWeight(ATOOLS::Flavour & flav)
    \brief If the flavour is found in the CCMap, its total weight (flavour weight times spin weight) is 
    returned; if not the method yields 0.
  */
  /*!
    \fn double Constituents::Smearing(ATOOLS::Flavour & flav)
    \brief If the flavour is found in the CCMap, its angular smearing parameter is returned; if not the method 
    yields 0.
  */
  /*!
    \fn double Constituents::ISpin(ATOOLS::Flavour & flav)
    \brief If the flavour is found in the CCMap, its integer spin is returned; if not the method yields 0.
  */
  /*!
    \fn double Constituents::PrintConstituents()
    \brief Prints the ingedients of the CCMap.
  */
}

#endif
