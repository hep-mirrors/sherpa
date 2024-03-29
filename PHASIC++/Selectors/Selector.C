#include "PHASIC++/Selectors/Selector.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Selector_Base
#define PARAMETER_TYPE PHASIC::Selector_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "PHASIC++/Process/Process_Base.H"

using namespace PHASIC;
using namespace ATOOLS;

void Selector_Log::Output() 
{ 
  msg_Info()<<"  Selector "<<m_name<<" rejection quota  : "
	    <<double(m_rejected)/double(m_rejected+m_passed)
	    <<"  ("<<m_rejected<<" / "<<m_passed+m_rejected<<")"<<std::endl;
}

std::vector<ATOOLS::Scoped_Settings> Selector_Key::GetSelectors() const
{
  Scoped_Settings addselsettings{ m_yaml };
  auto selectors = addselsettings["SELECTORS"].GetItems();
  auto userdefinedselectors = m_settings.GetItems();
  std::copy(selectors.begin(), selectors.end(),
            std::back_inserter(userdefinedselectors));
  return userdefinedselectors;
}

void Selector_Key::AddSelectorYAML(const std::string& yaml)
{
  if (m_yaml.empty())
    m_yaml = "SELECTORS:";
  m_yaml += "\n- " + yaml;
}

Selector_Base::Selector_Base(const std::string &name,Process_Base *const proc):
  m_name(name), m_on(false), m_isnlo(false),
  m_sel_log(new Selector_Log(m_name)), p_proc(proc),
  m_nin(p_proc?p_proc->NIn():0), m_nout(p_proc?p_proc->NOut():0),
  m_n(m_nin+m_nout), m_pass(1), p_sub(NULL),
  p_fl(p_proc?(Flavour*)&p_proc->Flavours().front():NULL),
  m_smin(0.), m_smax(sqr(rpa->gen.Ecms())),
  m_results{Weights_Map{1.0}}
{
  if (p_proc && p_proc->Info().Has(nlo_type::real|nlo_type::rsub))
    m_isnlo=true;
}

Selector_Base::~Selector_Base()
{
  if (m_sel_log!=NULL) delete m_sel_log;
}

bool Selector_Base::Trigger(const Vec4D_Vector &p,const Flavour *fl, size_t n)
{
  THROW(fatal_error,"Virtual function not reimplemented.");
  return false;
}

void Selector_Base::AddOnshellCondition(std::string,double)
{
}

void Selector_Base::Output() { 
  if (!(msg_LevelIsTracking())) return;
  if(m_sel_log) {
    m_sel_log->Output();
    msg_Out()<<m_name<<"  total number of rejections: "
	     <<m_sel_log->Rejections()<<std::endl;
  }
}

void Selector_Base::ReadInSubSelectors(const Selector_Key &key)
{
  for (auto s : key.m_settings[m_name]["Subselectors"].GetItems()) {
    Selector_Key subkey;
    subkey.m_settings = s;
    subkey.p_proc = key.p_proc;
    std::string type;
    if (s.IsList()) {
      type = s.SetDefault<std::string>({}).GetVector<std::string>()[0];
    } else {
      if (s.GetKeys().size() != 1)
        THROW(fatal_error,
              "Mapping-like selector settings must consist of exactly one"
              " key-value pair, where the key is the name of the selector,"
              " and the value is another mapping for the selector settings.");
      type = s.GetKeys().front();
    }
    auto* sel = Selector_Getter::GetObject(type, subkey);
    if (sel!=NULL) m_sels.push_back(sel);
    else THROW(fatal_error, "Did not find selector \""+type+"\".");
  }
}

const std::vector<ATOOLS::Weights_Map>& Selector_Base::Results() const
{
  return m_results;
}

void Selector_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  Selector_Getter::PrintGetterInfo(msg->Out(),25);
  msg_Out()<<"\n}"<<std::endl;
}

// default selector

namespace PHASIC {

  class No_Selector: public Selector_Base {
  public:

    No_Selector(): Selector_Base("No_Selector") {}

    bool Trigger(const Vec4D_Vector &,const Particle_List * pl=NULL) { return true; }
    bool Trigger(Selector_List &) { return true; }

    void BuildCuts(Cut_Data * cuts) {}

  };

}

DECLARE_ND_GETTER(No_Selector,"None",Selector_Base,Selector_Key,false);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,No_Selector>::
operator()(const Selector_Key &key) const
{
  return new No_Selector();
}

void ATOOLS::Getter<Selector_Base,Selector_Key,No_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"dummy selector"; 
}
