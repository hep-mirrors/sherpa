#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "ATOOLS/Org/Data_Reader.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Scale_Setter_Base
#define PARAMETER_TYPE PHASIC::Scale_Setter_Arguments
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

std::ostream &operator<<(std::ostream &ostr,const stp::id &scl)
{
  switch (scl) {
  case stp::ren: return ostr<<"ren";
  case stp::fac: return ostr<<"fac";
  case stp::size: return ostr<<"<error>";
  }
  return ostr<<"<unknown>";
}

Scale_Setter_Base::Scale_Setter_Base
(const Scale_Setter_Arguments &args,const bool scale2):
  p_proc(args.p_proc), p_model(args.p_model), p_cpls(args.p_cpls),
  m_scale(stp::size), m_coupling(args.m_coupling),
  m_scale2(scale2), p_jf(NULL)
{
}

void Scale_Setter_Base::SetCouplings()
{
  DEBUG_FUNC(p_proc->Name());
  if (p_cpls==NULL) THROW(fatal_error,"No coupling information");
  Data_Reader read(" ",",","#",":");
  std::vector<std::vector<std::string> > helpsvv;
  read.SetAddCommandLine(false);
  read.SetString(m_coupling);
  read.MatrixFromString(helpsvv,"");
  for (size_t i(0);i<helpsvv.size();++i) {
    if (helpsvv[i].size()!=2)
      THROW(fatal_error,"Invalid tag "+m_coupling+".");
    msg_Debugging()<<helpsvv[i][0]<<" -> "<<helpsvv[i][1]<<"\n";
    Coupling_Map::iterator cit(p_cpls->find(helpsvv[i][0]));
    if (cit!=p_cpls->end()) {
      if (ToType<size_t>(helpsvv[i][1])>=m_scale.size())
	THROW(fatal_error,"Index too large for "+helpsvv[i][0]+".");
      cit->second->SetScale(&m_scale[ToType<int>(helpsvv[i][1])]);
    }
    else {
      msg_Error()<<METHOD<<"("<<p_proc->Name()<<"): Valid tags are\n ";
      for (Coupling_Map::const_iterator cit(p_cpls->begin());
	   cit!=p_cpls->end();++cit) msg_Error()<<" "<<cit->first;
      msg_Error()<<"\n";
      THROW(fatal_error,"Invalid coupling tag "+helpsvv[i][0]+".");
    }
  }
}

Scale_Setter_Base::~Scale_Setter_Base()
{
}

double Scale_Setter_Base::CalculateScale2(const std::vector<ATOOLS::Vec4D> &p)
{
  return m_scale[stp::fac];
}

void Scale_Setter_Base::ShowSyntax(const size_t i)
{
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD<<"(): {\n\n"
	   <<"   // available scale choices\n\n";
  Scale_Getter_Function::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

Vec4D Scale_Setter_Base::Momentum(const size_t &i) const
{
  if (i>p_proc->NIn()+p_proc->NOut())
    THROW(fatal_error,"Momentum index too large");
  return p_proc->Integrator()->Momenta()[i];
}

double Scale_Setter_Base::HT() const
{
  double ht(0.0);
  const Vec4D_Vector &p(p_proc->Integrator()->Momenta());
  for (size_t i(p_proc->NIn());i<p.size();++i) ht+=p[i].PPerp();
  return ht;
}

double Scale_Setter_Base::YCut() const
{
  if (p_proc->LookUp() && p_proc->IsMapped())
    return p_proc->MapProc()->ScaleSetter()->YCut();
  if (p_jf!=NULL) return p_jf->Ycut();
  p_jf=(Jet_Finder*)p_proc->Selector()->GetSelector("Jetfinder");
  if (p_jf==NULL) THROW(critical_error,"Jet finder not found");
  return p_jf->Ycut();
}
