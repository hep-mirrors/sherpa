#include "AddOns/Analysis/Observables/HT.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Analysis_Key &key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const double min{ s["Min"].SetDefault(0.0).Get<double>() };
  const double max{ s["Max"].SetDefault(1.0).Get<double>() };
  const size_t bins{ s["Bins"].SetDefault(100).Get<size_t>() };
  const std::string scale{
    s["Scale"].SetDefault("Lin").Get<std::string>() };
  const std::string list{
    s["List"].SetDefault(std::string{finalstate_list}).Get<std::string>() };
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key &key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 0, Max: 1, Bins: 100, Scale: Lin, List: <list>}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_OBSERVABLE_GETTER(HT,"HT")
 
HT::HT(int type,double xmin,double xmax,int nbins,
       const std::string & listname,const std::string & reflistname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_reflist=reflistname;
  if (listname!="") {
    m_listname = listname;
    m_name = listname;
    if (m_reflist!="" && m_reflist!="FinalState") m_name += "_"+m_reflist;
    m_name+="_HT.dat";
  }
  else
    m_name = "HT.dat";
}

void HT::Evaluate(const ATOOLS::Particle_List& pl,
		  double weight, double ncount)
{
  ATOOLS::Particle_List* ref=p_ana->GetParticleList(m_reflist);
  ATOOLS::Particle_List* jets=p_ana->GetParticleList(m_listname);
  double HT=0.0;
  if(jets->size()==0 || ref==NULL || ref->empty()) {
    p_histo->Insert(0.0,0.0,ncount);
    return;
  }
  for (ATOOLS::Particle_List::const_iterator pit=jets->begin();
       pit!=jets->end();++pit) {
    HT+=(*pit)->Momentum().EPerp();
  }
  p_histo->Insert(HT,weight,ncount);
}

void HT::EvaluateNLOcontrib(double weight,double ncount )
{
  ATOOLS::Particle_List* ref=p_ana->GetParticleList(m_reflist);
  ATOOLS::Particle_List* jets=p_ana->GetParticleList(m_listname);
  double HT=0.0;
  if(jets->size()==0 || ref==NULL || ref->empty()) {
    p_histo->InsertMCB(0.0,0.0,ncount);
    return;
  }
  for (ATOOLS::Particle_List::const_iterator pit=jets->begin();
       pit!=jets->end();++pit) {
    HT+=(*pit)->Momentum().EPerp();
  }
  p_histo->InsertMCB(HT,weight,ncount);
}

void HT::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * HT::Copy() const 
{
  return new HT(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}
