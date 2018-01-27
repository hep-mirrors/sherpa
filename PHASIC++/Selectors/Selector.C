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

#include "PHASIC++/Process/Process_Base.H"

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
  DEBUG_FUNC("'"<<path<<"','"<<file<<"'");
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
  msg_Debugging()<<(*this);
}

std::ostream & PHASIC::operator<<(std::ostream & s,
                                  const PHASIC::Selector_Key & sk)
{
  s<<"Process: "<<(sk.p_proc?sk.p_proc->Name():"unknown")<<"\n";
  s<<"Key:     "<<sk.m_key<<", del="<<sk.m_del<<"\n";
  for (size_t i(0);i<sk.size();++i) s<<sk[i]<<std::endl;
  return s;
}

Selector_Base::Selector_Base(const std::string &name,Process_Base *const proc):
  m_name(name), m_on(false), m_isnlo(false),
  m_sel_log(new Selector_Log(m_name)), p_proc(proc),
  m_nin(p_proc?p_proc->NIn():0), m_nout(p_proc?p_proc->NOut():0),
  m_n(m_nin+m_nout), m_pass(1), p_sub(NULL),
  p_fl(p_proc?(Flavour*)&p_proc->Flavours().front():NULL),
  m_smin(0.), m_smax(sqr(rpa->gen.Ecms()))
{
  if (p_proc->Info().Has(nlo_type::real|nlo_type::rsub)) m_isnlo=true;
}

bool Selector_Base::RSTrigger(NLO_subevtlist *const subs)
{
  Flavour_Vector fl(p_fl,&p_fl[m_n]);
  int pass(0), nout(m_nout), nn(m_n);
  for (size_t n(0);n<subs->size();++n) {
    p_sub=(*subs)[n];
    m_nout=(m_n=p_sub->m_n)-m_nin;
    for (size_t i(0);i<m_n;++i) p_fl[i]=p_sub->p_fl[i];
    Vec4D_Vector mom(p_sub->p_mom,&p_sub->p_mom[m_n]);
    for (size_t i(0);i<m_nin;++i)
      if (mom[i][0]<0.0) mom[i]=-mom[i];
    Selector_List sl=Selector_List
      (p_sub->p_fl,p_sub->m_n,mom,m_nin);
    if (!Trigger(sl)) p_sub->m_trig=0;
    if (p_sub->m_trig) pass=1;
    p_sub=NULL;
  }
  m_n=nn;
  m_nout=nout;
  for (size_t i(0);i<m_n;++i) p_fl[i]=fl[i];
  return pass;
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

void Selector_Base::ReadInSubSelectors(const Selector_Key &key,size_t idx)
{
  DEBUG_FUNC("idx="<<idx);
  size_t open(0);
  Selector_Key * subkey(NULL);
  for (size_t k=idx;k<key.size();++k) {
    if (open==0) {
      if (subkey) THROW(fatal_error,"Read-in error.");
      subkey = new Selector_Key(key.p_proc,key.p_read);
      subkey->m_key=key[k][0];
    }
    if      (open==0 && key[k].back()=="{")                 open++;
    else if (open==1 && key[k].size()==1 && key[k][0]=="}") open--;
    else {
      if (open==0) {
        subkey->push_back(std::vector<std::string>(key[k].size()-1));
        for (size_t j=1;j<key[k].size();++j) subkey->back()[j-1]=key[k][j];
      }
      else {
        if      (key[k].back()=="{")                 open++;
        else if (key[k].size()==1 && key[k][0]=="}") open--;
        subkey->push_back(std::vector<std::string>(key[k].size()));
        for (size_t j=0;j<key[k].size();++j) subkey->back()[j]=key[k][j];
      }
    }
    if (open==0) {
      Selector_Base *sel(Selector_Getter::GetObject(subkey->m_key,*subkey));
      if (sel!=NULL) m_sels.push_back(sel);
      else THROW(fatal_error, "Did not find selector \""+subkey->m_key+"\".");
      if (msg_LevelIsDebugging()) {
        msg_Debugging()<<"subkey:\n";
        msg_Debugging()<<"  "<<subkey->m_key<<std::endl;
        for (size_t i(0);i<subkey->size();++i)
          msg_Debugging()<<"    "<<(*subkey)[i]<<std::endl;
        msg_Debugging()<<"-> found "<<(sel?sel->Name():"none")<<std::endl;
      }
      subkey=NULL;
    }
  }
  if (subkey) THROW(fatal_error,"Read-in error. Selector not processed.");
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
