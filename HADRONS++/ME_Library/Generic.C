#include "HADRONS++/ME_Library/Generic.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

#include "METOOLS/Main/Partial_Amplitude_Base.H"

Generic::Generic(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& decayindices,
                 const std::string& name):
  HD_ME_Base(flavs,decayindices,name)
{
  p_me=Partial_Amplitude_Base::Select(flavs);
  if(size()!=p_me->size())
    THROW(fatal_error, "size()!=p_me->size()");
}

Generic::~Generic() {
  delete p_me;
}

void Generic::Calculate(const Vec4D_Vector& p, bool m_anti)
{
  p_me->Calculate(p, m_anti);
  for(size_t i(0);i<size();++i) {
    (*this)[i]=m_factor*(*p_me)[i];
  }
}

DEFINE_ME_GETTER(Generic,"Generic")

void ATOOLS::Getter<HD_ME_Base,ME_Parameters,Generic>::
PrintInfo(std::ostream &str,const size_t width) const {
  str<<"Chooses a generic matrix element according to the spin structure.";
}
