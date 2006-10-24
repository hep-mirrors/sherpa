#include "Primitive_Observable_Base.H"
#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const GetOrdering(const String_Matrix &parameters)
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
  Primitive_Observable_Base *					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetOrdering<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"inlist outlist"; }

#define DEFINE_ORDERING_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;

namespace ANALYSIS {

  class Order_PT {
  public:
    bool operator()(const Particle *a,const Particle *b)
    {
      return a->Momentum().PPerp2()>b->Momentum().PPerp2();
    }
  };

  class PT_Ordering: public Primitive_Observable_Base {
  private:

    std::string m_outlist;

  public:

    PT_Ordering(const std::string &inlist,const std::string &outlist);
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight,int ncount);
    Primitive_Observable_Base *Copy() const;
    void EndEvaluation(double scale);

  };

  DEFINE_ORDERING_GETTER(PT_Ordering,PT_Ordering_Getter,"PTOrder");
    
  PT_Ordering::PT_Ordering(const std::string &inlist,
			   const std::string &outlist):
    m_outlist(outlist)
  {
    m_splitt_flag = false;
    m_listname=inlist;
  }

  void PT_Ordering::Evaluate(const ATOOLS::Particle_List &particlelist,
			     double weight,int ncount)
  {
    Particle_List *outlist(new Particle_List(particlelist.size()));
    p_ana->AddParticleList(m_outlist,outlist);
    for (size_t i(0);i<particlelist.size();++i) 
      (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
    std::sort(outlist->begin(),outlist->end(),Order_PT());
  }

  Primitive_Observable_Base *PT_Ordering::Copy() const
  {
    return new PT_Ordering(m_listname,m_outlist);
  }

  void PT_Ordering::EndEvaluation(double scale)
  {
  }
  
  class Order_ET {
  public:
    bool operator()(const Particle *a,const Particle *b)
    {
      return a->Momentum().EPerp()>b->Momentum().EPerp();
    }
  };

  class ET_Ordering: public Primitive_Observable_Base {
  private:

    std::string m_outlist;

  public:

    ET_Ordering(const std::string &inlist,const std::string &outlist);
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight,int ncount);
    Primitive_Observable_Base *Copy() const;
    void EndEvaluation(double scale);

  };

  DEFINE_ORDERING_GETTER(ET_Ordering,ET_Ordering_Getter,"ETOrder");
    
  ET_Ordering::ET_Ordering(const std::string &inlist,
			   const std::string &outlist):
    m_outlist(outlist)
  {
    m_splitt_flag = false;
    m_listname=inlist;
  }

  void ET_Ordering::Evaluate(const ATOOLS::Particle_List &particlelist,
			     double weight,int ncount)
  {
    Particle_List *outlist(new Particle_List());
    p_ana->AddParticleList(m_outlist,outlist);
    for (size_t i(0);i<particlelist.size();++i) 
      outlist->push_back(new ATOOLS::Particle(*particlelist[i]));
    std::sort(outlist->begin(),outlist->end(),Order_ET());
  }

  Primitive_Observable_Base *ET_Ordering::Copy() const
  {
    return new ET_Ordering(m_listname,m_outlist);
  }

  void ET_Ordering::EndEvaluation(double scale)
  {
  }
  
}
