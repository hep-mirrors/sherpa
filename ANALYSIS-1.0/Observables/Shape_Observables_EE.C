#include "Shape_Observables_EE.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"EEShapes";
    return new Class(10*(int)(parameters[0][3]=="Log"),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100, scale=0;
  std::string list="Analysed";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(scale,min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *const					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|Log [list] -> EEShapes"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME)

#include "Primitive_Analysis.H"

using namespace ATOOLS;
using namespace std;


Event_Shapes_Observable_Base::Event_Shapes_Observable_Base(int _type,double _min,double _max,int _nbins,
							   const string & _name) :
  Primitive_Observable_Base(_type,_min,_max,_nbins,NULL),
  m_key(std::string("EvtShapeData"))
{
  m_name = _name+string(".dat");
}


//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Thrust,Thrust_Getter,"Thrust");

Thrust::Thrust(int _type,double _min,double _max,int _nbins,const string & _name) :
  Event_Shapes_Observable_Base(_type,_min,_max,_nbins,_name) {}

void Thrust::Evaluate(const Particle_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(1.-data->Get<Event_Shape_EE_Data>().thrust,weight,ncount);
  }
}

Primitive_Observable_Base * Thrust::Copy() const 
{
  return new Thrust(m_type,m_xmin,m_xmax,m_nbins);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Major,Major_Getter,"Major");

Major::Major(int _type,double _min,double _max,int _nbins,const string & _name) :
  Event_Shapes_Observable_Base(_type,_min,_max,_nbins,_name) {}

void Major::Evaluate(const Particle_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().major,weight,ncount);
  }
}

Primitive_Observable_Base * Major::Copy() const 
{
  return new Major(m_type,m_xmin,m_xmax,m_nbins);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Minor,Minor_Getter,"Minor");

Minor::Minor(int _type,double _min,double _max,int _nbins,const string & _name) :
  Event_Shapes_Observable_Base(_type,_min,_max,_nbins,_name) {}

void Minor::Evaluate(const Particle_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().minor,weight,ncount);
  }
}

Primitive_Observable_Base * Minor::Copy() const 
{
  return new Minor(m_type,m_xmin,m_xmax,m_nbins);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Oblateness,Oblateness_Getter,"Oblat");

Oblateness::Oblateness(int _type,double _min,double _max,int _nbins,const string & _name) :
  Event_Shapes_Observable_Base(_type,_min,_max,_nbins,_name) {}

void Oblateness::Evaluate(const Particle_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().oblateness,weight,ncount);
  }
}

Primitive_Observable_Base * Oblateness::Copy() const 
{
  return new Oblateness(m_type,m_xmin,m_xmax,m_nbins);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(PT_In_Thrust,PT_In_Thrust_Getter,"PTIn");

PT_In_Thrust::PT_In_Thrust(int _type,double _min,double _max,int _nbins,
			   const std::string & _inlistname,const string & _name) :
  Event_Shapes_Observable_Base(_type,_min,_max,_nbins,_name) 
{ m_inlistname = _inlistname; }

void PT_In_Thrust::Evaluate(const Particle_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Majoraxis = data->Get<Event_Shape_EE_Data>().majoraxis;
    Particle_List * pl=p_ana->GetParticleList(m_inlistname);
    for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
      p_histo->Insert(dabs(Majoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

Primitive_Observable_Base * PT_In_Thrust::Copy() const 
{
  return new PT_In_Thrust(m_type,m_xmin,m_xmax,m_nbins,m_inlistname);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(PT_Out_Thrust,PT_Out_Thrust_Getter,"PTOut");

PT_Out_Thrust::PT_Out_Thrust(int _type,double _min,double _max,int _nbins,
			     const std::string & _inlistname,const string & _name) :
  Event_Shapes_Observable_Base(_type,_min,_max,_nbins,_name) 
{ m_inlistname = _inlistname; }


void PT_Out_Thrust::Evaluate(const Particle_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Minoraxis = data->Get<Event_Shape_EE_Data>().minoraxis;
    Particle_List * pl=p_ana->GetParticleList(m_inlistname);
    for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
      p_histo->Insert(dabs(Minoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

Primitive_Observable_Base * PT_Out_Thrust::Copy() const 
{
  return new PT_Out_Thrust(m_type,m_xmin,m_xmax,m_nbins,m_inlistname);
}

