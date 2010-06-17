#include "PHASIC++/Selectors/Combined_Selector.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Combined_Selector::Combined_Selector(Process_Integrator *const proc):
  Selector_Base("Combined_Selector"), m_count(0), m_update(false)
{
  p_proc=proc;
}

Combined_Selector::~Combined_Selector()
{
  while (m_sels.size()>0) {
    delete *m_sels.begin();
    m_sels.erase(m_sels.begin());
  }
}

Selector_Key Combined_Selector::FindArguments
(const Selector_Key &strings,size_t &starty,size_t &startx)
{
  size_t j=0, open=0, cdel=0;
  Selector_Key result(p_proc,strings.p_read);
  if (strings[starty].size()>startx) j=startx;
  for (size_t i=starty;i<strings.size();++i) {
    result.push_back(std::vector<std::string>(strings[i].size()-j));
    for (size_t k=0;j<strings[i].size();++j,++k) {
      result.back()[k]=strings[i][j];
      size_t opos=result.back()[k].find("{");
      if (opos!=std::string::npos) {
	++open;
	if (open==1 && k==0) {
	  cdel=1;
	  if (opos+1<result.back()[k].size()) 
	    result.back()[k]=result.back()[k].substr(opos+1);
	  else 
	    result.back()[k]="";
	  if (result.back()[k].length()==0) {
	    --k;
	    result.back().pop_back();
	    if (result.back().size()==0) result.pop_back();
	    continue;
	  }
	}
      }
      if (open>0) {
	size_t cpos=result.back()[k].find("}");
	if (cpos!=std::string::npos) --open;
	if (open==0 && cdel==1) {
	  result.back()[k]=result.back()[k].substr(0,cpos);
	  result.back().resize(k+1);
	  if (k==0 && result.back()[0].length()==0) 
	    result.resize(result.size()-1);
	  starty=i;
	  return result;
	}
      }
    }
    if (open==0) break;
    j=0;
  }  
  return result;
}

bool Combined_Selector::Initialize(const Selector_Key &key)
{
  msg_Debugging()<<METHOD<<"(): '"<<p_proc->Process()->Name()<<"' {\n";
  for (size_t k=0;k<key.size();++k) {
    if (key[k].size()>0) {
      if (key[k][0]=="{" || key[k][0]=="}") continue;
    }
    size_t col=1;
    size_t startk=k;
    Selector_Key mat(FindArguments(key,k,col));
    mat.m_key=key[startk][0];
    Selector_Base *sel(Selector_Getter::GetObject(key[startk][0],mat));
    if (sel!=NULL) {
      m_sels.push_back(sel);
      m_update|=sel->NeedUpdate();
      if (msg_LevelIsDebugging()) {
	msg_Out()<<"      new Selector_Base(\""
		 <<key[k][0]<<"\",";
	for (size_t i=0;i<mat.size();++i) {
	  msg_Out()<<"{"<<(mat[i].size()>0?mat[i][0]:"");
	  for (size_t j=1;j<mat[i].size();++j) 
	    msg_Out()<<","<<mat[i][j];
	  msg_Out()<<"}";
	}
	msg_Out()<<")\n";
      }
    }
    else {
      msg_Out()<<endl;
      THROW(fatal_error, "Did not find selector \""+key[startk][0]+"\".");
    }
  }
  msg_Debugging()<<"   }\n";
  return true;
}

int  Combined_Selector::NeedUpdate() 
{ 
  return m_update>0; 
}

bool Combined_Selector::Trigger(const Vec4D_Vector &p) 
{
  m_count++;
  for (size_t i=0; i<m_sels.size(); ++i) {
    if (!(m_sels[i]->Trigger(p))) return 0;
  }
  return 1;
}

bool Combined_Selector::JetTrigger
(const Vec4D_Vector &p,const Flavour_Vector &fl, int ns)
{
  for (size_t i=0; i<m_sels.size(); ++i) {
    if (!(m_sels[i]->JetTrigger(p,fl,ns))) return 0;
  }
  return 1;
}

bool Combined_Selector::JetTrigger(const Vec4D_Vector &p)
{
  for (size_t i=0; i<m_sels.size(); ++i) {
    if (!(m_sels[i]->JetTrigger(p))) return 0;
  }
  return 1;
}

bool Combined_Selector::NoJetTrigger(const Vec4D_Vector &p)
{
  m_count++;
  for (size_t i=0; i<m_sels.size(); ++i) {
    if (!(m_sels[i]->NoJetTrigger(p))) return 0;
  }
  return 1;
}

void Combined_Selector::BuildCuts(Cut_Data * cuts)
{
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (!m_sels[i]->IsConditional()) m_sels[i]->BuildCuts(cuts);

  //smin update!!!

  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->IsConditional()) m_sels[i]->BuildCuts(cuts);
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->IsConditional()) m_sels[i]->BuildCuts(cuts);

  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
  cuts->Complete();
  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
}

void Combined_Selector::UpdateCuts(double sprime,double y,Cut_Data * cuts)
{
  cuts->Reset(NeedUpdate());
  if (NeedUpdate()) 
    for (size_t i=0; i<m_sels.size(); ++i) m_sels[i]->UpdateCuts(sprime,y,cuts);
  cuts->Update(sprime,y);      
  for (size_t i=0; i<m_osc.size(); ++i) cuts->Setscut(m_osc[i].first,m_osc[i].second);
}

 
void Combined_Selector::AddOnshellCondition(std::string s,double d)
{
  m_osc.push_back(std::pair<std::string,double>(s,d));
}

void Combined_Selector::Output()
{
  msg_Debugging()<<"========================================="<<std::endl
			    <<"Efficiency of the Selector : "<<m_name<<std::endl;
  for (size_t i=0; i<m_sels.size(); ++i) m_sels[i]->Output();
  msg_Debugging()<<"========================================="<<std::endl;
}

Selector_Base * Combined_Selector::GetSelector(const std::string &name) const
{
  for (size_t i=0; i<m_sels.size(); ++i) 
    if (m_sels[i]->Name()==name) return m_sels[i];
  return 0;
  
}

