#include "AddOns/Analysis/Observables/ChargedDEta.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() < 6)
    THROW(missing_input, "Missing parameter values.");
  const auto list = parameters.size() > 6 ? parameters[6] : finalstate_list;
  return new Class(HistogramType(parameters[5]),
                   s.Interprete<double>(parameters[2]),
                   s.Interprete<double>(parameters[3]),
                   s.Interprete<int>(parameters[4]),list,
                   s.Interprete<int>(parameters[0]),
                   s.Interprete<int>(parameters[1]));
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"[kf1, kf2, min, max, bins, Lin|LinErr|Log|LogErr, list] ... list is optional"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_OBSERVABLE_GETTER(ChargedDEta,ChargedDEta_Getter,"ChargedDEta")
 
ChargedDEta::ChargedDEta(int type,double xmin,double xmax,int nbins,
                         const std::string & listname,
                         const int kf1, const int kf2) :
Primitive_Observable_Base(type,xmin,xmax,nbins), m_flav1(kf1), m_flav2(kf2)
{
  m_name = "ChargedDEta_"+m_flav1.ShellName()+m_flav2.ShellName()+".dat";
  m_listname = listname;
}

void ChargedDEta::Evaluate(const ATOOLS::Particle_List& pl,
                           double weight, double ncount)
{
  ATOOLS::Particle_List* list=p_ana->GetParticleList(m_listname);
  std::vector<ATOOLS::Vec4D> p1;
  std::vector<double> charges1;
  std::vector<ATOOLS::Vec4D> p2;
  for (ATOOLS::Particle_List::const_iterator it=list->begin();
       it!=list->end(); ++it) {
    if (m_flav1.Includes((*it)->Flav())) {
      p1.push_back((*it)->Momentum());
      charges1.push_back((*it)->Flav().Charge());
    }
    else if (m_flav2.Includes((*it)->Flav())) p2.push_back((*it)->Momentum());
  }
  for (size_t i=0; i<p1.size(); ++i) {
    for (size_t j=0; j<p2.size(); ++j) {
      double DEta=p2[j].DEta(p1[i]);
      p_histo->Insert(charges1[i]*DEta,weight,ncount);
    }
  }
}


Primitive_Observable_Base * ChargedDEta::Copy() const 
{
  return new ChargedDEta(m_type,m_xmin,m_xmax,m_nbins,m_listname,int(m_flav1),
                         int(m_flav2));
}
