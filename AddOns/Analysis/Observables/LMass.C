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
Primitive_Observable_Base *GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(0.0).Get<double>();
  const auto max = s["Max"].SetDefault(1.0).Get<double>();
  const auto bins = s["Bins"].SetDefault(100).Get<size_t>();
  const auto scale = s["Scale"].SetDefault("Lin").Get<std::string>();
  const auto list = s["List"]
    .SetDefault(std::string(finalstate_list))
    .Get<std::string>();
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 1, Max: 10, Bins: 100, Scale: Lin, List: FinalState}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
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
