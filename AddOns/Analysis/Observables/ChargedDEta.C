#include "AddOns/Analysis/Observables/ChargedDEta.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<6) return NULL;
    std::string list=parameters[0].size()>6?parameters[0][6]:finalstate_list;
    return new Class(HistogramType(parameters[0][5]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),list,
                     ATOOLS::ToType<int>(parameters[0][0]),
                     ATOOLS::ToType<int>(parameters[0][1]));
  }
  else {
    return NULL;
  }
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,CLASS>::operator()(const Argument_Matrix &parameters) const \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"kf1 kf2 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Argument_Matrix);	\
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
