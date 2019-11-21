#include "AddOns/Analysis/Observables/CD_Parameters.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Message.H"

namespace ANALYSIS {
  std::ostream& operator<<( std::ostream& ostr, const CD_Parameter_Data & data) {
    ostr<<"CD_Parameter_Data : "<<data.cparameter<<","<<data.dparameter;
    return ostr;
  }
}

using namespace ANALYSIS;
using namespace ATOOLS;

#include "ATOOLS/Org/MyStrStream.H"

DECLARE_GETTER(CD_Parameter_Calculator,"CDCalc",
	       Analysis_Object,Analysis_Key);

Analysis_Object * ATOOLS::Getter<Analysis_Object,Analysis_Key,
				 CD_Parameter_Calculator>::operator()(const Analysis_Key& key) const
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto listname = s.SetDefault(finalstate_list).Get<std::string>();
  return new CD_Parameter_Calculator(listname);
}

void ATOOLS::Getter<Analysis_Object,Analysis_Key,CD_Parameter_Calculator>::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list";
}

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

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 0, Max: 1, Bins: 100, Scale: Lin, List: FinalState} ... depends on CDCalc"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)


CD_Parameter_Calculator::CD_Parameter_Calculator(const std::string & listname)
  : m_tensor(1), m_key(listname+"_CD_Parameters") 
{
  m_name = listname+"CD_Parameters_Calculator";
  m_listname = listname;
}

void CD_Parameter_Calculator::Evaluate(const Blob_List & ,double weight, double ncount) {
  Particle_List * pl = p_ana->GetParticleList(m_listname);
  if (pl==NULL) {
    msg_Out()<<"WARNING in CD_Parameter_Calculator::Evaluate : particle list "<<m_listname<<" not found "<<std::endl;
    return;
  }
  double lambda1=0., lambda2=0., lambda3=0.;
  if (pl->size()>0) {
    m_tensor.Calculate(*pl);
    lambda1=m_tensor.EigenValue(0);
    lambda2=m_tensor.EigenValue(1);
    lambda3=m_tensor.EigenValue(2);
  }
  double cparameter= 3.*(lambda1*lambda2+lambda2*lambda3+lambda3*lambda1);
  double dparameter= 27.*lambda1*lambda2*lambda3;

  p_ana->AddData(m_key,
		 new Blob_Data<CD_Parameter_Data>(CD_Parameter_Data(cparameter,dparameter)));
}

void CD_Parameter_Calculator::EvaluateNLOcontrib(double value, double ncount)
{
  Blob_List bl;
  Evaluate(bl,value, ncount);
}

void CD_Parameter_Calculator::EvaluateNLOevt()
{ }

Analysis_Object * CD_Parameter_Calculator::GetCopy() const
{
  return new CD_Parameter_Calculator(m_listname);
}



// ----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(C_Parameter,
			 C_Parameter_Getter,"CParam")


C_Parameter::C_Parameter(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_CD_Parameters")
{
  m_listname = listname;
  m_name = std::string("C_Parameter.dat");
}

void C_Parameter::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<CD_Parameter_Data>().cparameter,weight,ncount);
  }
}
 
void C_Parameter::EvaluateNLOcontrib(double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<CD_Parameter_Data>().cparameter,weight,ncount);
  }
}
 
void C_Parameter::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}
 
Primitive_Observable_Base * C_Parameter::Copy() const
{
  return new C_Parameter(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}



// ----------------------------------------------------------------------


DEFINE_OBSERVABLE_GETTER(D_Parameter,
			 D_Parameter_Getter,"DParam")


D_Parameter::D_Parameter(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_CD_Parameters")
{
  m_listname = listname;
  m_name = std::string("D_Parameter.dat");
}

void D_Parameter::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<CD_Parameter_Data>().dparameter,weight,ncount);
  }
}
 
void D_Parameter::EvaluateNLOcontrib(double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<CD_Parameter_Data>().dparameter,weight,ncount);
  }
}
 
void D_Parameter::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}
 
Primitive_Observable_Base * D_Parameter::Copy() const
{
  return new D_Parameter(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

namespace ATOOLS {

template <> Blob_Data<CD_Parameter_Data>::~Blob_Data() { }
template class Blob_Data<CD_Parameter_Data>;

}
