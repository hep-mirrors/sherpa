#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include <map>

namespace ANALYSIS {

  struct kfitem {
    long int kf;
    size_t item;
  };

  class Universal_Observable: public Primitive_Observable_Base, public ATOOLS::Tag_Replacer {  
  private:
    ATOOLS::Algebra_Interpreter m_interpreter;
    std::string  m_expression,m_oexp;
    bool m_init;
    std::map<std::string,kfitem> m_itemlist;

    bool CheckInList();
  public:

    Universal_Observable(const std::string &obs,int type,double xmin,double xmax,int nbins,
       const std::string & listname=std::string(""));
    
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    

  };// end of class Universal_Observable

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace std;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>5?parameters[0][5]:finalstate_list;
    return new Class(parameters[0][0],HistogramType(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  else if (parameters.size()<5) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, obs="0.", scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
    else if (parameters[i][0]=="OBSERVABLE") obs=parameters[i][1];
  }
  return new Class(obs,HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"observable min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

DEFINE_OBSERVABLE_GETTER(Universal_Observable,Universal_Observable_Getter,"UniObs")
 
Universal_Observable::Universal_Observable(const std::string &obs,int type,double xmin,double xmax,
					   int nbins,const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_oexp=m_expression=obs;
  if (listname!="") {
    m_listname = listname;
    m_name = "UniObs_"+listname+".dat";
  }
  else
    m_name = "UniObs.dat";

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
  m_interpreter.AddTag("P_SUM","(10.,1.,4.,5.)");
  pos=m_expression.find("item[-");
  while (pos!=string::npos) {
    m_expression[pos+5]='m';
    pos=m_expression.find("item[-");
  }
}

bool Universal_Observable::CheckInList()
{
  ATOOLS::Particle_List * plist=p_ana->GetParticleList(m_listname);
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

void Universal_Observable::Evaluate(const ATOOLS::Particle_List& pl,
				    double weight, double ncount)
{
  if (!CheckInList()) {
    p_histo->Insert(0.,0.,ncount);
    return;
  }

  if (!m_init) {
    m_interpreter.Interprete(m_expression);
    m_init=1;
  }

  double value=(m_interpreter.Calculate())->Get<double>();
  p_histo->Insert(value,weight,ncount);
}

void Universal_Observable::EvaluateNLOcontrib(double weight,double ncount )
{
  if (!CheckInList()) {
    p_histo->InsertMCB(0.,0.,ncount);
    return;
  }

  if (!m_init) {
    m_interpreter.Interprete(m_expression);
    m_init=1;
  }

  double value=(m_interpreter.Calculate())->Get<double>();

  p_histo->InsertMCB(value,weight,ncount);
}

void Universal_Observable::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Universal_Observable::Copy() const 
{
  return new Universal_Observable(m_oexp,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

std::string Universal_Observable::ReplaceTags(std::string &expr) const
{
  return m_interpreter.ReplaceTags(expr);
}

ATOOLS::Term *Universal_Observable::ReplaceTags(ATOOLS::Term *term) const
{
  ATOOLS::Particle_List * plist=p_ana->GetParticleList(m_listname);

  if (term->Tag()=="H_T") {
    double ht(0.0);
    for  (size_t i=0;i<plist->size();++i)
      ht+=(*plist)[i]->Momentum().PPerp();
    term->Set(ht);
    return term;
  }

  if (term->Tag()=="P_SUM") {
    ATOOLS::Vec4D ps(0.,0.,0.,0.);
    for  (size_t i=0;i<plist->size();++i)
      ps+=(*plist)[i]->Momentum();
    term->Set(ps);
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
