#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>
#include <map>

namespace ANALYSIS {

  struct kfitem {
    long int kf;
    size_t item;
  };

  class Uni_Selector: public Trigger_Base, public ATOOLS::Tag_Replacer {  
  protected:
    ATOOLS::Algebra_Interpreter m_interpreter;
    std::string  m_expression;
    double m_xmin, m_xmax;
    bool m_init;
    std::map<std::string,kfitem> m_itemlist;

    bool CheckInList();
  public:

    Uni_Selector
    (const std::string &expr,
     const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    Analysis_Object *GetCopy() const;
  };// end of class Uni_Selector


}// end of namespace ANALYSIS

using namespace ANALYSIS;
using namespace std;

template <class Class>
Analysis_Object *const 
GetUniSelector(const Argument_Matrix &parameters) 
{				
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    return new Class(parameters[0][0],
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     parameters[0][3],parameters[0][4]);
  }
  if (parameters.size()<5) return NULL;
  double min=30.0, max=70.0; 
  std::string expr="", inlist="", outlist="";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Expression") expr=parameters[i][1];
  }
  return new Class(expr,min,max,inlist,outlist);
}									

#define DEFINE_UNI_SELECTOR_GETTER_METHOD(CLASS,NAME)	\
  Analysis_Object *				\
  NAME::operator()(const Argument_Matrix &parameters) const	\
  { return GetUniSelector<CLASS>(parameters); }

#define DEFINE_UNI_SELECTOR_PRINT_METHOD(NAME)			\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"expression min max inlist outlist"; }

#define DEFINE_UNI_SELECTOR_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_UNI_SELECTOR_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_UNI_SELECTOR_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_UNI_SELECTOR_GETTER(Uni_Selector,
				 Uni_Selector_Getter,"UniSel")


Uni_Selector::
Uni_Selector(const std::string &expr,
		   const double min,const double max,
		   const std::string &inlist,
		   const std::string &outlist):
  Trigger_Base(inlist,outlist)
{
  m_expression=expr;
  m_xmin=min;
  m_xmax=max;
  m_init=0;
  m_interpreter.SetTagReplacer(this);
  string help=m_expression;
  size_t pos=help.find("item[");
  while (pos!=string::npos) {
    help=help.substr(pos);
    size_t psep=help.find(":");
    size_t pend=help.find("]");
    string sitem=help.substr(0,pend+1);
    if (sitem[5]=='-') sitem[5]='m';
    if (m_itemlist.find(sitem)==m_itemlist.end()) {
      kfitem hki;
      hki.kf=kf_none;
       if (psep!=string::npos&&psep<pend) {
	hki.kf=ATOOLS::ToType<int>(help.substr(5,psep-5));
      }
      else psep=4;
      hki.item=ATOOLS::ToType<size_t>(help.substr(psep+1,pend-psep-1));
      m_itemlist[sitem]=hki;
      m_interpreter.AddTag(sitem,"(5.,0.,3.,4.)");
    }
    help=help.substr(pend+1);
    pos=help.find("item[");
  }
  m_interpreter.AddTag("H_T","1.0");
  pos=m_expression.find("item[-");
  while (pos!=string::npos) {
    m_expression[pos+5]='m';
    pos=m_expression.find("item[-");
  }
}

bool Uni_Selector::CheckInList()
{
  ATOOLS::Particle_List * plist=p_ana->GetParticleList(m_inlist);
  if (plist->size()==0) return 0;

  for (map<string,kfitem>::iterator mit=m_itemlist.begin();mit!=m_itemlist.end();mit++) {
    kfitem hkfi=mit->second;
    
    int no=-1; 
    size_t pos=std::string::npos;
    for (size_t i=0;i<plist->size();++i) {
      if ((long int)((*plist)[i]->Flav())==hkfi.kf || hkfi.kf==kf_none) {
	++no;
	if (no==(int)hkfi.item) {
	  pos=i;
	  break;
	}
      }
    }
    if (pos==std::string::npos) return 0;
  }
  return 1;
}

void Uni_Selector::Evaluate(const ATOOLS::Particle_List &inlist,
				  ATOOLS::Particle_List &outlist,
				  double weight,double ncount)
{
  if (!CheckInList()) return;

  if (!m_init) {
    m_interpreter.Interprete(m_expression);
    m_init=1;
  }
  double value=(m_interpreter.Calculate())->Get<double>();
  if (value<m_xmin||value>m_xmax) return;

  outlist.resize(inlist.size());
  for (size_t i=0;i<inlist.size();++i) 
    outlist[i] = new ATOOLS::Particle(*inlist[i]);
}

Analysis_Object *Uni_Selector::GetCopy() const
{
  return new Uni_Selector(m_expression,m_xmin,m_xmax,
				m_inlist,m_outlist);
}

std::string Uni_Selector::ReplaceTags(std::string &expr) const
{
  return m_interpreter.ReplaceTags(expr);
}

ATOOLS::Term *Uni_Selector::ReplaceTags(ATOOLS::Term *term) const
{
  ATOOLS::Particle_List * plist=p_ana->GetParticleList(m_inlist);

  if (term->Tag()=="H_T") {
    double ht(0.0);
    for  (size_t i=0;i<plist->size();++i)
      ht+=(*plist)[i]->Momentum().PPerp();
    term->Set(ht);
    return term;
  }

  kfitem hkfi=m_itemlist.find(term->Tag())->second;
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<plist->size();++i) {
    if ((long int)((*plist)[i]->Flav())==hkfi.kf || hkfi.kf==kf_none) {
      ++no;
      if (no==(int)hkfi.item) {
	pos=i;
	break;
      }
    }
  }
  if (pos==std::string::npos) term->Set(ATOOLS::Vec4D(0.,0.,0.,0.));
  else term->Set((*plist)[pos]->Momentum());
  return term;
}
