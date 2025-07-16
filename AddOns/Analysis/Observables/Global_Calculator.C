#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Total_Momentum: public Primitive_Observable_Base {  
  protected:

    std::string  m_outlist;

  public:

    Total_Momentum(const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Total_Momentum

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *
GetGlobalCalculator(const Analysis_Key& key)
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size()<2) return NULL;
  return new Class(parameters[0],parameters[1]);
  return NULL;
}

#define DEFINE_GLOBAL_CALCULATOR_GETTER_METHOD(CLASS)		\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key &key) const \
  { return GetGlobalCalculator<CLASS>(key); }

#define DEFINE_GLOBAL_CALCULATOR_PRINT_METHOD(NAME)		\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"[inlist, resulttag]"; }

#define DEFINE_GLOBAL_CALCULATOR_GETTER(CLASS,TAG)		\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GLOBAL_CALCULATOR_GETTER_METHOD(CLASS)		\
  DEFINE_GLOBAL_CALCULATOR_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_GLOBAL_CALCULATOR_GETTER(Total_Momentum,"MomSum")

Total_Momentum::Total_Momentum(const std::string &inlist,
			       const std::string &outlist):
  m_outlist(outlist)
{
  m_splitt_flag = false;
  m_listname=inlist;
}

void Total_Momentum::Evaluate(const ATOOLS::Particle_List &inlist,
			      double weight,double ncount)
{
  Vec4D sum;
  for (size_t i(0);i<inlist.size();++i) sum+=inlist[i]->Momentum();
  msg_Debugging()<<"adding '"<<m_outlist<<"' = "<<sum<<"\n";
  p_ana->AddData(m_outlist,new Blob_Data<Vec4D>(sum));
}

Primitive_Observable_Base *Total_Momentum::Copy() const
{
  return new Total_Momentum(m_listname,m_outlist);
}

