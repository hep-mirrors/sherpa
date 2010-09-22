#include "PHASIC++/Selectors/Selector.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Selector_Base
#define PARAMETER_TYPE PHASIC::Selector_Key
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;

void Selector_Log::Output() 
{ 
  msg_Info()<<"  Selector "<<m_name<<" rejection quota  : "
	    <<double(m_rejected)/double(m_rejected+m_passed)
	    <<"  ("<<m_rejected<<" / "<<m_passed+m_rejected<<")"<<std::endl;
}

Selector_Key::~Selector_Key()
{
  if (m_del) delete p_read;
}

void Selector_Key::SetData(const std::string &tag,
			   const std::vector<std::string> &args)
{
  bool found(false);
  for (std::vector<std::vector<std::string> >::iterator 
	 lit(begin());lit!=end();++lit) {
    if (lit->front()==tag) {
      if (found) lit=erase(lit);
      else { 
	found=true;
	lit->resize(1);
	lit->insert(lit->end(),args.begin(),args.end());
      }
    }
  }
  if (!found) {
    std::vector<std::string> line(1,tag);
    line.insert(line.end(),args.begin(),args.end());
    push_back(line);
  }
}

void Selector_Key::ReadData(const std::string &path,const std::string &file)
{
  msg_Debugging()<<METHOD<<"('"<<path<<"','"<<file<<"'): {\n";
  msg_Indent();
  if (m_del && p_read!=NULL) delete p_read;
  p_read=new Data_Reader(" ",";","!");
  p_read->AddWordSeparator("\t");
  p_read->AddComment("#");
  p_read->AddComment("//");
  p_read->SetAddCommandLine(false);
  p_read->SetInputPath(path);
  p_read->SetInputFile(file);
  p_read->SetMatrixType(mtc::transposed);
  p_read->MatrixFromFile((*this),"");
  msg_Debugging()<<"}"<<std::endl;
}


Selector_Base::~Selector_Base() 
{ 
  if (m_sel_log!=NULL) delete m_sel_log;
}

bool Selector_Base::JetTrigger
(const Vec4D_Vector &,NLO_subevtlist *const sub)
{
  THROW(fatal_error,"Virtual method not redefined");
  return false;
}

bool Selector_Base::NoJetTrigger(const Vec4D_Vector &)
{
  THROW(fatal_error,"Virtual method not redefined");
  return false;
}

void Selector_Base::BuildCuts(Cut_Data *) 
{
  THROW(fatal_error,"Virtual method not redefined");
}

void Selector_Base::AddOnshellCondition(std::string,double) 
{
}

int Selector_Base::IsConditional()
{ 
  return 0; 
}

void Selector_Base::Output() { 
  if (!(msg_LevelIsTracking())) return;
  if(m_sel_log) {
    m_sel_log->Output();
    msg_Out()<<m_name<<"  total number of rejections: "
	     <<m_sel_log->Rejections()<<std::endl;
  }
}

void Selector_Base::ShowSyntax(const int mode)
{
  if (!msg_LevelIsInfo() || mode==0) return;
  msg_Out()<<METHOD<<"(): {\n\n";
  Selector_Getter::PrintGetterInfo(msg->Out(),20);
  msg_Out()<<"\n}"<<std::endl;
}

// default selector

namespace PHASIC {

  class No_Selector: public Selector_Base {
  public:

    No_Selector(): Selector_Base("No_Selector") {}

    bool Trigger(const Vec4D_Vector &) { return true; }
    bool JetTrigger(const Vec4D_Vector &,NLO_subevtlist *const) 
    { return true; }
    bool NoJetTrigger(const Vec4D_Vector &) { return true; }

    void SetRange(ATOOLS::Flavour_Vector,double,double) {}
    void BuildCuts(Cut_Data * cuts) {}
    void Output() {}

  };

}

DECLARE_ND_GETTER(No_Selector_Getter,"None",Selector_Base,Selector_Key,false);

Selector_Base *No_Selector_Getter::operator()(const Selector_Key &key) const
{
  return new No_Selector();
}

void No_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"dummy selector"; 
}
