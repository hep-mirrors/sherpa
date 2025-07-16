#include "AddOns/Analysis/Observables/Sphericity_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

namespace ANALYSIS {
  std::ostream& operator<<( std::ostream& ostr, const Sphericity_Data & data) {
    ostr<<"Sphericity_Data : "<<data.sphericity<<","<<data.aplanarity<<","<<data.planarity;
    return ostr;
  }
}

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Sphericity_Calculator,"SphCalc",
	       Analysis_Object,Analysis_Key);

Analysis_Object * ATOOLS::Getter<Analysis_Object,Analysis_Key,
				 Sphericity_Calculator>::operator()(const Analysis_Key& key) const
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto listname = s.SetDefault(finalstate_list).Get<std::string>();
  return new Sphericity_Calculator(listname);
}

void ATOOLS::Getter<Analysis_Object,Analysis_Key,Sphericity_Calculator>::PrintInfo(std::ostream &str,const size_t width) const
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
  { str<<"e.g. {Min: 0, Max: 1, Bins: 100, Scale: Lin, List: FinalState} ... depends on SphCalc"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)


Sphericity_Calculator::Sphericity_Calculator(const std::string & listname)
  : m_tensor(2), m_key(listname+"_Sphericity_Data") 
{
  m_name = std::string("Sphericitys_Calculator");
  m_listname = listname;
}

void Sphericity_Calculator::Evaluate(const Blob_List & ,double weight, double ncount) {
  Particle_List * pl = p_ana->GetParticleList(m_listname);
  if (pl==NULL) {
    msg_Out()<<"WARNING in Sphericity_Calculator::Evaluate : particle list "<<m_listname<<" not found "<<std::endl;
    return;
  }
  double lambda2=0., lambda3=0.;
  if (pl->size()>0) {
    m_tensor.Calculate(*pl);
    lambda2=m_tensor.EigenValue(1);
    lambda3=m_tensor.EigenValue(2);
  }
  double sphericity = 3./2.*(lambda2+lambda3);
  double aplanarity = 3./2.*lambda3;
  double planarity  = lambda2-lambda3;

  p_ana->AddData(m_key,
		 new Blob_Data<Sphericity_Data>(Sphericity_Data(sphericity,aplanarity,planarity)));
}

void Sphericity_Calculator::EvaluateNLOcontrib(double value, double ncount)
{
  Blob_List bl;
  Evaluate(bl,value, ncount);
}

void Sphericity_Calculator::EvaluateNLOevt()
{ }

Analysis_Object * Sphericity_Calculator::GetCopy() const
{
  return new Sphericity_Calculator(m_listname);
}



// ----------------------------------------------------------------------


DEFINE_OBSERVABLE_GETTER(Sphericity,
			 Sphericity_Getter,"Sphericity")


Sphericity::Sphericity(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_Sphericity_Data")
{
  m_listname = listname;
  m_name = std::string("Sphericity.dat");
}

void Sphericity::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Sphericity_Data>().sphericity,weight,ncount);
  }
}
 
 
void Sphericity::EvaluateNLOcontrib(double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<Sphericity_Data>().sphericity,weight,ncount);
  }
}
 
void Sphericity::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Sphericity::Copy() const
{
  return new Sphericity(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}



// ----------------------------------------------------------------------


DEFINE_OBSERVABLE_GETTER(Aplanarity,
			 Aplanarity_Getter,"Aplanarity")


Aplanarity::Aplanarity(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_Sphericity_Data")
{
  m_listname = listname;
  m_name = std::string("Aplanarity.dat");
}

void Aplanarity::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Sphericity_Data>().aplanarity,weight,ncount);
  }
}
 
void Aplanarity::EvaluateNLOcontrib(double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<Sphericity_Data>().aplanarity,weight,ncount);
  }
}
 
void Aplanarity::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Aplanarity::Copy() const
{
  return new Aplanarity(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


// ----------------------------------------------------------------------


DEFINE_OBSERVABLE_GETTER(Planarity,
			 Planarity_Getter,"Planarity")


Planarity::Planarity(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_Sphericity_Data")
{
  m_listname = listname;
  m_name = std::string("Planarity.dat");
}

void Planarity::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Sphericity_Data>().planarity,weight,ncount);
  }
}
 
 
void Planarity::EvaluateNLOcontrib(double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<Sphericity_Data>().planarity,weight,ncount);
  }
}
 
void Planarity::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Planarity::Copy() const
{
  return new Planarity(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

namespace ATOOLS {

template <> Blob_Data<Sphericity_Data>::~Blob_Data() { }
template class Blob_Data<Sphericity_Data>;

}
