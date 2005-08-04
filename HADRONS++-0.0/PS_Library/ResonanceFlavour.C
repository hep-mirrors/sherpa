#include "Flavour.H"
#include "ResonanceFlavour.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

ResonanceFlavour::ResonanceFlavour( string _name, double _mass, double _width )
{
  name = _name;
  mass = _mass;
  width = _width;
}
ResonanceFlavour::ResonanceFlavour( kf::code _kfc, double _mass, double _width )
{
  name = Flavour(_kfc).Name();
  mass = _mass;
  width = _width;
}
 
void ResonanceFlavour::Set( kf::code _kfc, double _mass, double _width )
{
  name = Flavour(_kfc).Name();
  mass = _mass;
  width = _width;
}
