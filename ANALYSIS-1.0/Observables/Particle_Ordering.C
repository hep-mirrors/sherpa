#include "Primitive_Observable_Base.H"
#include "Primitive_Analysis.H"
#include "MyStrStream.H"
#include <algorithm>

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const GetOrdering(const Argument_Matrix &parameters)
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
  Primitive_Observable_Base *						\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetOrdering<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"inlist outlist"; }

#define DEFINE_ORDERING_GETTER(CLASS,NAME,TAG)				\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
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
  class NAME: public Primitive_Observable_Base {			\
  private:								\
    std::string m_outlist;						\
  public:								\
    NAME(const std::string &inlist,const std::string &outlist);		\
    void Evaluate(const ATOOLS::Particle_List &particlelist,		\
		  double weight,int ncount);				\
    Primitive_Observable_Base *Copy() const;				\
    void EndEvaluation(double scale);					\
  };									\
  DEFINE_ORDERING_GETTER(NAME,GNAME,TAG)				\
  NAME::NAME(const std::string &inlist,const std::string &outlist):	\
    m_outlist(outlist)							\
  {									\
    m_splitt_flag = false;						\
    m_listname=inlist;							\
  }									\
  void NAME::Evaluate(const ATOOLS::Particle_List &particlelist,	\
		      double weight,int ncount)				\
  {									\
    Particle_List *outlist(new Particle_List());			\
    outlist->resize(particlelist.size());				\
    for (size_t i(0);i<particlelist.size();++i)				\
      (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);		\
    std::sort(outlist->begin(),outlist->end(),CNAME());			\
    p_ana->AddParticleList(m_outlist,outlist);				\
  }									\
  Primitive_Observable_Base *NAME::Copy() const				\
  {									\
    return new NAME(m_listname,m_outlist);				\
  }									\
  void NAME::EndEvaluation(double scale)				\
  {									\
  }									
  

namespace ANALYSIS {

  DEFINE_PARTICLE_ORDERING(Order_PT,Order_PT_Getter,Sort_PT,PPerp2,"PTOrder")
  DEFINE_PARTICLE_ORDERING(Order_ET,Order_ET_Getter,Sort_ET,EPerp,"ETOrder")
  DEFINE_PARTICLE_ORDERING(Order_Y,Order_Y_Getter,Sort_Y,Y,"YOrder")
  DEFINE_PARTICLE_ORDERING(Order_Eta,Order_Eta_Getter,Sort_Eta,Eta,"EtaOrder")
  DEFINE_PARTICLE_ORDERING(Order_Phi,Order_Phi_Getter,Sort_Phi,Phi,"PhiOrder")
  
}
