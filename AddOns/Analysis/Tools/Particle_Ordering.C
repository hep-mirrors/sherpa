#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;

template <class Class>
Analysis_Object *const GetOrdering(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<2) return NULL;
    return new Class(parameters[0][0],parameters[0][1]);
  }
  else if (parameters.size()<2) return NULL;
  std::string inlist="FinalState", outlist="Ordered";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="INLIST") inlist=parameters[i][1];
    else if (parameters[i][0]=="OUTLIST") outlist=parameters[i][1];
  }
  return new Class(inlist,outlist);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Analysis_Object *							\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetOrdering<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"inlist outlist"; }

#define DEFINE_ORDERING_GETTER(CLASS,NAME,TAG)				\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);		\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;

#define DEFINE_PARTICLE_ORDERING(NAME,GNAME,CNAME,VAR,TAG)		\
  class CNAME {								\
  public:								\
    bool operator()(const Particle *a,const Particle *b)		\
    {									\
      return a->Momentum().VAR()>b->Momentum().VAR();			\
    }									\
  };									\
  class NAME: public Analysis_Object {					\
  private:								\
    std::string m_inlist, m_outlist;					\
  public:								\
    NAME(const std::string &inlist,const std::string &outlist);		\
    void Evaluate(const ATOOLS::Blob_List &bl,				\
		  double weight,double ncount);				\
    Analysis_Object *GetCopy() const;					\
  };									\
  DEFINE_ORDERING_GETTER(NAME,GNAME,TAG)				\
  NAME::NAME(const std::string &inlist,const std::string &outlist):	\
    m_inlist(inlist), m_outlist(outlist) {}				\
  void NAME::Evaluate(const ATOOLS::Blob_List &bl,			\
		      double weight,double ncount)				\
  {									\
    Particle_List *outlist(new Particle_List());			\
    ATOOLS::Particle_List *inlist(p_ana->GetParticleList(m_inlist));	\
    if (inlist==NULL) {							\
      msg_Error()<<METHOD<<"(): List '"<<m_inlist		\
			 <<"' not found."<<std::endl;			\
      p_ana->AddParticleList(m_outlist,outlist);			\
      return;								\
    }									\
    outlist->resize(inlist->size());					\
    for (size_t i(0);i<inlist->size();++i)				\
      (*outlist)[i] = new ATOOLS::Particle(*(*inlist)[i]);		\
    std::sort(outlist->begin(),outlist->end(),CNAME());			\
    p_ana->AddParticleList(m_outlist,outlist);				\
  }									\
  Analysis_Object *NAME::GetCopy() const				\
  {									\
    return new NAME(m_inlist,m_outlist);				\
  }								       
  

namespace ANALYSIS {

  DEFINE_PARTICLE_ORDERING(Order_PT,Order_PT_Getter,Sort_PT,PPerp2,"PTOrder")
  DEFINE_PARTICLE_ORDERING(Order_ET,Order_ET_Getter,Sort_ET,EPerp,"ETOrder")
  DEFINE_PARTICLE_ORDERING(Order_Y,Order_Y_Getter,Sort_Y,Y,"YOrder")
  DEFINE_PARTICLE_ORDERING(Order_Eta,Order_Eta_Getter,Sort_Eta,Eta,"EtaOrder")
  DEFINE_PARTICLE_ORDERING(Order_Phi,Order_Phi_Getter,Sort_Phi,Phi,"PhiOrder")
  
}
