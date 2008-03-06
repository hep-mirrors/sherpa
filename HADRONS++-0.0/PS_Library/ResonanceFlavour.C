#include "Flavour.H"
#include "ResonanceFlavour.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

ResonanceFlavour::ResonanceFlavour( string _name, double _mass, double _width )
{
  m_name = _name;
  m_mass = _mass;
  m_width = _width;
}
ResonanceFlavour::ResonanceFlavour( kf_code _kfc, double _mass, double _width )
{
  m_name = Flavour(_kfc).IDName();
  m_mass = _mass;
  m_width = _width;
}
 
void ResonanceFlavour::Set( kf_code _kfc, double _mass, double _width )
{
  m_name = Flavour(_kfc).IDName();
  m_mass = _mass;
  m_width = _width;
}
