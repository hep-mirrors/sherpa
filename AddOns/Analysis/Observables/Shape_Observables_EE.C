#include "AddOns/Analysis/Observables/Shape_Observables_EE.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

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
  const auto list = s["List"].SetDefault("EEShapes").Get<std::string>();
  return new Class(HistogramType(scale),min,max,bins,list);
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,CLASS>::operator()(const Analysis_Key& key) const \
  { return GetObservable<CLASS>(key); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Analysis_Key,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 1, Max: 10, Bins: 100, Scale: Lin, List: EEShapes} ... depends on EEShapes"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Analysis_Key);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ATOOLS;
using namespace std;


Event_Shapes_Observable_Base::Event_Shapes_Observable_Base(int _type,double _min,double _max,int _nbins,
							   const string & _name) :
  Primitive_Observable_Base(_type,_min,_max,_nbins),
  m_key(std::string("EvtShapeData"))
{
  m_name = _name+string(".dat");
}


//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Thrust,Thrust_Getter,"Thrust")

Thrust::Thrust(int type, double xmin, double xmax, int nbins, 
	       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{
}

void Thrust::Evaluate(const Blob_List &,double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(1.-data->Get<Event_Shape_EE_Data>().thrust,weight,ncount);
  }
}

void Thrust::EvaluateNLOcontrib(double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(1.-data->Get<Event_Shape_EE_Data>().thrust,weight,ncount);
  }
}

void Thrust::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Thrust::Copy() const 
{
  return new Thrust(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Major,Major_Getter,"Major")

Major::Major(int type, double xmin, double xmax, int nbins, 
		       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{
}

void Major::Evaluate(const Blob_List &,double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().major,weight,ncount);
  }
}

void Major::EvaluateNLOcontrib(double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<Event_Shape_EE_Data>().major,weight,ncount);
  }
}

void Major::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Major::Copy() const 
{
  return new Major(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Minor,Minor_Getter,"Minor")

Minor::Minor(int type, double xmin, double xmax, int nbins, 
		       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{
}

void Minor::Evaluate(const Blob_List &,double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().minor,weight,ncount);
  }
}

void Minor::EvaluateNLOcontrib(double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<Event_Shape_EE_Data>().minor,weight,ncount);
  }
}

void Minor::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Minor::Copy() const 
{
  return new Minor(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Oblateness,Oblateness_Getter,"Oblat")

Oblateness::Oblateness(int type, double xmin, double xmax, int nbins, 
		       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{
}

void Oblateness::Evaluate(const Blob_List &,double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().oblateness,weight,ncount);
  }
}

void Oblateness::EvaluateNLOcontrib(double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->InsertMCB(data->Get<Event_Shape_EE_Data>().oblateness,weight,ncount);
  }
}

void Oblateness::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * Oblateness::Copy() const 
{
  return new Oblateness(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(PT_In_Thrust,PT_In_Thrust_Getter,"PTIn")

PT_In_Thrust::PT_In_Thrust(int type, double xmin, double xmax, int nbins, 
		       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{ 
  m_listname = listname; 
}

void PT_In_Thrust::Evaluate(const Particle_List & pl,double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Majoraxis = data->Get<Event_Shape_EE_Data>().majoraxis;
    for (Particle_List::const_iterator pit=pl.begin();pit!=pl.end();++pit) {
      p_histo->Insert(dabs(Majoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

void PT_In_Thrust::EvaluateNLOcontrib(double weight,double ncount) 
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Majoraxis = data->Get<Event_Shape_EE_Data>().majoraxis;
    for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
      p_histo->InsertMCB(dabs(Majoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

void PT_In_Thrust::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * PT_In_Thrust::Copy() const 
{
  return new PT_In_Thrust(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(PT_Out_Thrust,PT_Out_Thrust_Getter,"PTOut")

PT_Out_Thrust::PT_Out_Thrust(int type, double xmin, double xmax, int nbins, 
		       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{
  m_listname = listname; 
}


void PT_Out_Thrust::Evaluate(const Particle_List & pl,double weight,double ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Minoraxis = data->Get<Event_Shape_EE_Data>().minoraxis;
    for (Particle_List::const_iterator pit=pl.begin();pit!=pl.end();++pit) {
      p_histo->Insert(dabs(Minoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

void PT_Out_Thrust::EvaluateNLOcontrib(double weight,double ncount) 
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Minoraxis = data->Get<Event_Shape_EE_Data>().minoraxis;
    for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
      p_histo->InsertMCB(dabs(Minoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

void PT_Out_Thrust::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}

Primitive_Observable_Base * PT_Out_Thrust::Copy() const 
{
  return new PT_Out_Thrust(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Eta_Thrust,Eta_Thrust_Getter,"EtaThrust")

Eta_Thrust::Eta_Thrust(int type, double xmin, double xmax, int nbins, 
		       const std::string & listname, const std::string & name) :
  Event_Shapes_Observable_Base(type,xmin,xmax,nbins,name) 
{ 
  m_listname = listname; 
}


void Eta_Thrust::Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D thrust = data->Get<Event_Shape_EE_Data>().thrustaxis;
    double eta = Vec4D(0.,thrust).Eta();
    if (0<=m_xmin) eta=dabs(eta);
    p_histo->Insert(eta,weight,ncount);
  }
}

Primitive_Observable_Base * Eta_Thrust::Copy() const 
{
  return new Eta_Thrust(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

