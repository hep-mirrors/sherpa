#include "HADRONS++/ME_Library/Exchange_ME_Shim.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;
using namespace HADRONS;
using namespace std;

Exchange_ME_Shim::Exchange_ME_Shim(const Flavour_Vector& flavs,
                                   const vector<int>& indices,
                                   const string& name):
  HD_ME_Base(flavs,indices,name), p_me(NULL)
{
  // 1: a decay. Leg 0 is the decaying particle and everything else is
  // outgoing. This is Exchange_ME's default, but stated explicitly because
  // it is what makes the q sign and the topology checks correct here.
  p_me = new METOOLS::Exchange_ME(flavs,indices,name,1);
}

Exchange_ME_Shim::~Exchange_ME_Shim()
{
  if (p_me) delete p_me;
}

void Exchange_ME_Shim::SetModelParameters(Scoped_Settings& s)
{
  p_me->SetModelParameters(s);
  if (p_me->size()!=size())
    THROW(fatal_error,"Exchange_ME and its shim disagree on the number of "
          "spin combinations for "+Name()+". Both derive from "
          "Spin_Amplitudes over the same flavours, so this means the two "
          "were built with different index lists.");
}

void Exchange_ME_Shim::Calculate(const Vec4D_Vector& moms, bool anti)
{
  p_me->Calculate(moms,anti);
  // Both sides are Spin_Amplitudes over the same flavours, so the helicity
  // ordering is shared and this is a straight copy, not a remapping.
  for (size_t i(0);i<size();++i) (*this)[i] = (*p_me)[i];
}

// NOT via DEFINE_ME_GETTER. Both METOOLS/ME_Base.H and
// HADRONS++/ME_Library/HD_ME_Base.H define a macro of that name, against
// DIFFERENT registries, and this file includes both -- so whichever header
// lands last silently decides which registry the class is registered in.
// Spelled out, it cannot be decided by include order.
DECLARE_GETTER(HADRONS::Exchange_ME_Shim,"Exchange_ME",
               HADRONS::HD_ME_Base,METOOLS::ME_Parameters);

HADRONS::HD_ME_Base * ATOOLS::Getter<HADRONS::HD_ME_Base,
                                     METOOLS::ME_Parameters,
                                     HADRONS::Exchange_ME_Shim>::
operator()(const METOOLS::ME_Parameters &parameters) const
{ return new HADRONS::Exchange_ME_Shim(parameters.flavs,parameters.indices,
                                       "Exchange_ME"); }

void ATOOLS::Getter<HADRONS::HD_ME_Base,METOOLS::ME_Parameters,
                    HADRONS::Exchange_ME_Shim>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Exchange_ME: sum over exchanged bosons joining two currents.\n\n"
    <<"Declared as a Currents: block with an Exchange: list; an entry with\n"
    <<"no Boson: is a contact term, which is what Current_ME always was.\n";
}
