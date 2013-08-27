#include "PHOTONS++/Main/Dipole_Type.H"

using namespace std;
using namespace ATOOLS;
using namespace PHOTONS;

std::ostream& PHOTONS::operator<<(std::ostream& s,const IdPair& id)
{
  return s<<"("<<id.i<<","<<id.j<<")";
}

std::ostream& PHOTONS::operator<<(std::ostream& s,const IdPairNbar& idn)
{
  return s<<"["<<idn.ij<<": "<<idn.nbar<<"]";
}

