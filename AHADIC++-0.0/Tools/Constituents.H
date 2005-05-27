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
}

#endif
