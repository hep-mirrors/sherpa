#include "Shape_Observables_EE.H"
#include "Primitive_Analysis.H"

using namespace ANALYSIS;
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

