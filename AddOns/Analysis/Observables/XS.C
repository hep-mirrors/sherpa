#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {

  class XS: public Primitive_Observable_Base {  
  private:

  public:
    XS(const std::string & listname=std::string(""));
    
    void EvaluateNLOcontrib(double weigXS, double ncount);
    void EvaluateNLOevt();
    void Evaluate(const ATOOLS::Particle_List & pl, double weigXS, double ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class XS

}// end of namespace ANALYSIS

using namespace ANALYSIS;


template <class Class>
Primitive_Observable_Base* GetObservable(const Analysis_Key &key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto list = s["List"].SetDefault("FinalState").Get<std::string>();
  return new Class(list);
}


#define DEFINE_GETTER_METHOD(CLASS)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key &key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {List: FinalState}"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS)					\
  DEFINE_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_OBSERVABLE_GETTER(XS,"XS")

XS::XS(const std::string & listname) :
  Primitive_Observable_Base(101,0.,1.,1)
{
  if (listname!="") {
    m_listname = listname;
    m_name = listname;
    m_name+="_XS.dat";
  }
  else
    m_name = "XS.dat";
}

void XS::Evaluate(const ATOOLS::Particle_List& pl,
                  double weigXS, double ncount)
{
  p_histo->Insert(0.5,weigXS,ncount);
}

void XS::EvaluateNLOcontrib(double weigXS,double ncount )
{
  p_histo->InsertMCB(0.5,weigXS,ncount);
}

void XS::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * XS::Copy() const 
{
  return new XS(m_listname);
}
