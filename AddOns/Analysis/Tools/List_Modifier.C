#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Two_Particle_Modifier_Base: public Analysis_Object,
				    public Tag_Replacer {  
  protected:

    std::string m_crit, m_inlist, m_outlist;

    ATOOLS::Flavour m_flav[2];
    size_t          m_item[2];

    ATOOLS::Algebra_Interpreter m_interpreter;

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;        

  public:

    Two_Particle_Modifier_Base
    (const ATOOLS::Flavour flav[2],const size_t item[2],
     const std::string &crit,
     const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Blob_List &blobs,
		  double value=1.0,double ncount=1);
    
    virtual bool Modify(Particle *&p1,Particle *&p2) const = 0;

  };// end of class Two_Particle_Modifier_Base

  class Two_Particle_Swap: public Two_Particle_Modifier_Base {  
  public:

    Two_Particle_Swap(const ATOOLS::Flavour flav[2],const size_t item[2],
		      const std::string &crit,
		      const std::string &inlist,const std::string &outlist);
    
    bool Modify(Particle *&p1,Particle *&p2) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class Two_Particle_Swap

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetTwoParticleModifier(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    size_t item[2];
    ATOOLS::Flavour flav[2];
    for (size_t i(0);i<2;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][2*i]);
      flav[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) flav[i]=flav[i].Bar();
      item[i]=ATOOLS::ToType<size_t>(parameters[0][2*i+1]);
    }
    return new Class(flav,item,
		     parameters[0][4],parameters[0][5],parameters[0][6]);
  }
  return NULL;
}									

#define DEFINE_TWO_MODIFIER_GETTER_METHOD(CLASS,NAME)			\
  Analysis_Object *							\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetTwoParticleModifier<CLASS>(parameters); }

#define DEFINE_TWO_MODIFIER_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav1 item1 flav2 item2 crit inlist outlist"; }

#define DEFINE_TWO_MODIFIER_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_TWO_MODIFIER_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_TWO_MODIFIER_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Two_Particle_Modifier_Base::
Two_Particle_Modifier_Base(const ATOOLS::Flavour flav[2],const size_t item[2],
			   const std::string &crit,const std::string &inlist,
			   const std::string &outlist):
  m_crit(crit),
  m_outlist(outlist)
{
  for (size_t i(0);i<2;++i) {
    m_flav[i]=flav[i];
    m_item[i]=item[i];
  }
  m_inlist=inlist;
  msg_Debugging()<<METHOD<<"(): m_crit = '"<<m_crit<<"' {\n";
  size_t bpos(m_crit.find('@')); 
  while (bpos!=std::string::npos) {
    size_t epos(m_crit.find('@',bpos+1));
    if (epos!=std::string::npos && epos-bpos>1) {
      if (m_crit[epos-1]==']') {
	m_interpreter.AddTag(m_crit.substr(bpos,epos-bpos+1),"(1.0,0.0,0.0,1.0)");
	msg_Debugging()<<"  adding vector tag '"
		       <<m_crit.substr(bpos+1,epos-bpos-1)<<"'\n";
      }
      else {
	m_interpreter.AddTag(m_crit.substr(bpos,epos-bpos+1),"1");
	msg_Debugging()<<"  adding double tag '"
		       <<m_crit.substr(bpos+1,epos-bpos-1)<<"'\n";
      }
    }
    bpos=m_crit.find('@',++epos);
  }
  m_interpreter.SetTagReplacer(this);
  m_interpreter.Interprete(m_crit);
  msg_Debugging()<<"}\n";
}

void Two_Particle_Modifier_Base::Evaluate
(const ATOOLS::Blob_List &blobs,double value,double ncount)
{
  ATOOLS::Particle_List *inlist=p_ana->GetParticleList(m_inlist);
  if (inlist==NULL) {
    msg_Error()<<METHOD<<"(): List '"<<m_inlist
		       <<"' not found."<<std::endl;
    return;
  }
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  int no(-1);
  size_t pos[2]={std::string::npos,std::string::npos};
  for (size_t k(0);k<2;++k) {
    no=-1;
    for (size_t i(0);i<inlist->size();++i) {
      if ((*inlist)[i]->Flav()==m_flav[k] || 
	  m_flav[k].Kfcode()==kf_none) {
	++no;
	if (no==(int)m_item[k]) {
	  pos[k]=i;
	  break;
	}
      }
    }
  }
  if (pos[0]==std::string::npos || pos[1]==std::string::npos) return;
  if (m_interpreter.Calculate()->Get<double>()!=0.0) {
    outlist->resize(inlist->size());
    for (size_t i=0;i<inlist->size();++i) 
      (*outlist)[i] = new ATOOLS::Particle(*(*inlist)[i]);
    Modify((*outlist)[pos[0]],(*outlist)[pos[1]]);
  }
}

std::string Two_Particle_Modifier_Base::ReplaceTags(std::string &expr) const
{
  return m_interpreter.ReplaceTags(expr);
}

Term *Two_Particle_Modifier_Base::ReplaceTags(Term *term) const
{
  std::string tag(term->Tag().substr(1,term->Tag().length()-2));
  Blob_Data_Base *data((*p_ana)[tag]);
  if (data==NULL) THROW(critical_error,"Data '"+tag+"' not found.");
  if (term->Tag()[term->Tag().length()-1]==']') 
    term->Set(data->Get<Vec4D>());
  else term->Set(data->Get<double>());
  return term;
}

DEFINE_TWO_MODIFIER_GETTER(Two_Particle_Swap,
			   Two_Particle_Swap_Getter,"TwoSwap")
  
Two_Particle_Swap::
Two_Particle_Swap(const ATOOLS::Flavour flav[2],const size_t item[2],
		  const std::string &crit,
		  const std::string &inlist,const std::string &outlist):
  Two_Particle_Modifier_Base(flav,item,crit,inlist,outlist) {}

bool Two_Particle_Swap::Modify(Particle *&p1,Particle *&p2) const
{
  std::swap<Particle*>(p2,p1);
  return true;
}

Analysis_Object *Two_Particle_Swap::GetCopy() const
{
  return new Two_Particle_Swap(m_flav,m_item,m_crit,m_inlist,m_outlist);
}

