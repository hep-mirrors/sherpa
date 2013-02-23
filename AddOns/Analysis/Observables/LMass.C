#include "AddOns/Analysis/Main/Primitive_Analysis.H"

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class ListMass: public Primitive_Observable_Base {  
  public:

    ListMass(int type,double xmin,double xmax,int nbins,
       const std::string & listname=std::string(""));
    
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class ListMass

}// end of namespace ANALYSIS

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:finalstate_list;
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,CLASS>::operator()(const Argument_Matrix &parameters) const \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_OBSERVABLE_GETTER(ListMass,"ListMass")
 
ListMass::ListMass(int type,double xmin,double xmax,int nbins,
       const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  if (listname!="") {
    m_listname = listname;
    m_name = listname;
    m_name+="_ListMass.dat";
  }
  else
    m_name = "ListMass.dat";
}

void ListMass::Evaluate(const ATOOLS::Particle_List& pl,
		  double weight, double ncount)
{
  ATOOLS::Particle_List* jets=p_ana->GetParticleList(m_listname);
  ATOOLS::Vec4D momsum(0.,0.,0.,0.);
  for (ATOOLS::Particle_List::const_iterator pit=jets->begin();
       pit!=jets->end();++pit) {
    momsum+=(*pit)->Momentum();
  }
  p_histo->Insert(momsum.Mass(),weight,ncount);
}

void ListMass::EvaluateNLOcontrib(double weight,double ncount )
{
  ATOOLS::Particle_List* jets=p_ana->GetParticleList(m_listname);
  ATOOLS::Vec4D momsum(0.,0.,0.,0.);
  for (ATOOLS::Particle_List::const_iterator pit=jets->begin();
       pit!=jets->end();++pit) {
    momsum+=(*pit)->Momentum();
  }
  p_histo->InsertMCB(momsum.Mass(),weight,ncount);
}

void ListMass::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * ListMass::Copy() const 
{
  return new ListMass(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
