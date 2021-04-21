#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Tree_ME2_Base
#define PARAMETER_TYPE PHASIC::External_ME_Args
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;

Tree_ME2_Base::Tree_ME2_Base(const External_ME_Args& args) :
  p_aqcd(NULL), p_aqed(NULL), m_namps(0), m_norm(1.0), 
  m_flavs(args.Flavours())
{
}

Tree_ME2_Base::~Tree_ME2_Base() {}

std::vector<Complex> Tree_ME2_Base::GetAmplitudes(const size_t &id)
{
  return std::vector<Complex>();
}

Complex Tree_ME2_Base::GetPhase(const size_t &id)
{
  return Complex(1.0,0.0);
}

Complex Tree_ME2_Base::GetHelicityPhase(const Vec4D &pijt,const Vec4D &eps1)
{
  THROW(not_implemented,"Missing phase for interference term");
  return Complex(0.0,0.0);
}

std::vector<Tree_ME2_Base::Map_Info>
Tree_ME2_Base::GetFlavourHelicityMap()
{
  return std::vector<Map_Info>();
}

void Tree_ME2_Base::FillCombinations
(std::set<std::pair<size_t,size_t> > &combs,
 std::map<size_t,ATOOLS::Flavour_Vector> &fls)
{
}

void Tree_ME2_Base::CompleteCombinations
(std::set<std::pair<size_t,size_t> > &combs,
 std::map<size_t,ATOOLS::Flavour_Vector> &fls) const
{
  std::set<std::pair<size_t,size_t> > tcombs(combs);
  for (std::set<std::pair<size_t,size_t> >::const_iterator
	 cit=tcombs.begin();cit!=tcombs.end();++cit) {
    size_t ida(cit->first), idb(cit->second);
    size_t idc((1<<m_flavs.size())-1-ida-idb);
    msg_Debugging()<<"comb "<<ID(ida)
		   <<" "<<ID(idb)<<" "<<ID(idc)<<"\n";
    combs.insert(std::pair<size_t,size_t>(idb,ida));
    combs.insert(std::pair<size_t,size_t>(idb,idc));
    combs.insert(std::pair<size_t,size_t>(idc,idb));
    combs.insert(std::pair<size_t,size_t>(idc,ida));
    combs.insert(std::pair<size_t,size_t>(ida,idc));
  }
  std::map<size_t,ATOOLS::Flavour_Vector> tfls(fls);
  for (std::map<size_t,ATOOLS::Flavour_Vector>::const_iterator
	 cit=tfls.begin();cit!=tfls.end();++cit) {
    msg_Debugging()<<"flav "<<ID(cit->first)
		   <<" -> "<<cit->second<<"\n";
    Flavour_Vector fl(cit->second);
    for (size_t i(0);i<fl.size();++i) fl[i]=fl[i].Bar();
    fls[(1<<m_flavs.size())-1-cit->first]=fl;
  }
}

int Tree_ME2_Base::OrderQCD(const int &id) const
{
  THROW(not_implemented, "Not implemented");
}

int Tree_ME2_Base::OrderEW(const int &id) const
{
  THROW(not_implemented, "Not implemented");
}

double Tree_ME2_Base::TR() const
{
  return 0.5;
}

typedef ATOOLS::Getter_Function
<Tree_ME2_Base,PHASIC::External_ME_Args> Tree_ME2_Getter;

Tree_ME2_Base* Tree_ME2_Base::GetME2(const PHASIC::External_ME_Args& args)
{
  Tree_ME2_Getter::Getter_List glist(Tree_ME2_Getter::GetGetters());
  for (Tree_ME2_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    Tree_ME2_Base *me2=(*git)->GetObject(args);
    if (me2) return me2;
  }
  return NULL;
}

Tree_ME2_Base *Tree_ME2_Base::GetME2(const std::string& tag,
				     const External_ME_Args& args)
{
  Tree_ME2_Base* me2=Tree_ME2_Getter::GetObject(tag, args);
  if (me2==NULL) {
    THROW(fatal_error, "Did not find ME^2 "+tag);
  }
  else return me2;
}

Tree_ME2_Base* Tree_ME2_Base::GetME2(const PHASIC::Process_Info& pi)
 {
   External_ME_Args args(pi.m_ii.GetExternal(),
			 pi.m_fi.GetExternal(),
			 pi.m_maxcpl,pi.m_megenerator);
   return GetME2(args);
 }

void Tree_ME2_Base::SetCouplings(const MODEL::Coupling_Map& cpls)
{
  p_aqcd=cpls.Get("Alpha_QCD");
  p_aqed=cpls.Get("Alpha_QED");
}

double Tree_ME2_Base::AlphaQCD() const
{
  return p_aqcd ? p_aqcd->Default()*p_aqcd->Factor() : 
    s_model->ScalarConstant("alpha_S");
}

double Tree_ME2_Base::AlphaQED() const
{
  return p_aqed ? p_aqed->Default()*p_aqed->Factor() : 
    s_model->ScalarConstant("alpha_QED");
}

double Trivial_Tree::Calc(const ATOOLS::Vec4D_Vector &p)
{
  return 0.0;
}

namespace PHASIC {

  std::ostream &operator<<(std::ostream &str,const Tree_ME2_Base::Map_Info &mi)
  {
    return str<<'{'<<mi.m_id<<"|P"<<mi.m_perm<<",H"<<mi.m_hels<<'}';
  }

}
