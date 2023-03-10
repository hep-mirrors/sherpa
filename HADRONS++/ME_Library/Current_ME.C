#include "HADRONS++/ME_Library/Current_ME.H"
#include "HADRONS++/Current_Library/Current_Base.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

Current_ME::Current_ME(const ATOOLS::Flavour_Vector& flavs,
                       const std::vector<int>& decayindices,
                       const std::string& name):
  HD_ME_Base(flavs,decayindices,name), p_c1(NULL), p_c2(NULL)
{
}

void Current_ME::SetModelParameters(ATOOLS::Scoped_Settings& s)
{
  DEBUG_FUNC("");
  p_c1 = SelectCurrent(s["J1"]);
  p_c2 = SelectCurrent(s["J2"]);

  m_factor=Tools::GF/sqrt(2.);

  if(m_flavs.size() != p_c1->DecayIndices().size()+p_c2->DecayIndices().size())
    THROW(fatal_error, "Current selection does not look sane for "+Name());
}

Current_Base* Current_ME::SelectCurrent(ATOOLS::Scoped_Settings s)
{
  std::string current_string = s["Type"].SetDefault("Unknown").Get<string>();
  auto indices = s["Indices"].SetDefault({-1}).GetVector<int>();
  ME_Parameters fi(m_flavs, indices);
  Current_Base* current=Current_Getter_Function::GetObject(current_string,fi);
  if(current==NULL) {
    THROW(fatal_error, "Current '"+current_string+"' not found.");
  }

  GeneralModel model;
  for (const auto& key: s.GetKeys()) {
    if (key=="Type") continue;
    if (key=="Indices") continue;
    if (s[key].GetItemsCount()>1) {
      model[key+string("_abs")] = s[key][0].GetScalarWithOtherDefault(-1.0);
      model[key+string("_phase")] = s[key][1].GetScalarWithOtherDefault(0.0);
    }
    else model[key] = s[key].GetScalarWithOtherDefault(-1.0);
  }
  current->SetModelParameters(model);
  return current;
}


Current_ME::~Current_ME() {
  if (p_c1) delete p_c1;
  if (p_c2) delete p_c2;
}

void Current_ME::Calculate(const Vec4D_Vector& p, bool anti)
{
  p_c1->Calc(p, anti);
  p_c2->Calc(p, anti);
  
  std::vector<int> spins,spins1,spins2;
  for(size_t i=0;i<size();i++) {
    spins=GetSpinCombination(i);
    spins1.clear(); spins2.clear();
    for(size_t j=0;j<p_c1->DecayIndices().size();j++)
      spins1.push_back(spins[p_c1->DecayIndices()[j]]);
    for(size_t j=0;j<p_c2->DecayIndices().size();j++)
      spins2.push_back(spins[p_c2->DecayIndices()[j]]);
    // now we know the spin combinations in both currents
    // let's fill the results:
    (*this)[i]=m_factor*p_c1->Get(spins1)*p_c2->Get(spins2);
  }
}

DEFINE_ME_GETTER(Current_ME,"Current_ME")

void ATOOLS::Getter<HD_ME_Base,ME_Parameters,Current_ME>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Current_ME"<<endl;
}
