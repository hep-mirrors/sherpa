#include "PHOTONS++/Tools/Weight_Base.H"

using namespace PHOTONS;
using namespace std;

Weight_Base::Weight_Base() {
  m_weight    = 1;
  m_maxweight = 1;
  m_dtype     = Dipole_Type::unknown;
}

Weight_Base::~Weight_Base() {

}
