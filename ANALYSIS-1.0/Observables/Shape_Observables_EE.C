#include "Shape_Observables_EE.H"
#include "MyStrStream.H"


using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  return new Class(parameters);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list] -> EEShapes"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

#include "Primitive_Analysis.H"

using namespace ATOOLS;
using namespace std;


Event_Shapes_Observable_Base::Event_Shapes_Observable_Base(const String_Matrix & parameters) :
  Primitive_Observable_Base(parameters),
  m_key(std::string("EvtShapeData"))
{ }


//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Thrust,Thrust_Getter,"Thrust")

Thrust::Thrust(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{   
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Thrust in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Thrust.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Thrust";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "1-Thrust";
      p_histo->SetTitle(title);
    }
  }
}

Thrust::Thrust(const Thrust * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void Thrust::Evaluate(const Blob_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(1.-data->Get<Event_Shape_EE_Data>().thrust,weight,ncount);
  }
}

Primitive_Observable_Base * Thrust::Copy() const 
{
  return new Thrust(this);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Major,Major_Getter,"Major")

Major::Major(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{ 
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

Major::Major(const Major * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void Major::Evaluate(const Blob_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().major,weight,ncount);
  }
}

Primitive_Observable_Base * Major::Copy() const 
{
  return new Major(this);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Minor,Minor_Getter,"Minor")

Minor::Minor(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{ 
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

Minor::Minor(const Minor * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void Minor::Evaluate(const Blob_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().minor,weight,ncount);
  }
}

Primitive_Observable_Base * Minor::Copy() const 
{
  return new Minor(this);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Oblateness,Oblateness_Getter,"Oblat")

Oblateness::Oblateness(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{ 
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

Oblateness::Oblateness(const Oblateness * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void Oblateness::Evaluate(const Blob_List &,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<Event_Shape_EE_Data>().oblateness,weight,ncount);
  }
}

Primitive_Observable_Base * Oblateness::Copy() const 
{
  return new Oblateness(this);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(PT_In_Thrust,PT_In_Thrust_Getter,"PTIn")

PT_In_Thrust::PT_In_Thrust(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{ 
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

PT_In_Thrust::PT_In_Thrust(const PT_In_Thrust * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void PT_In_Thrust::Evaluate(const Particle_List & pl,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Majoraxis = data->Get<Event_Shape_EE_Data>().majoraxis;
    for (Particle_List::const_iterator pit=pl.begin();pit!=pl.end();++pit) {
      p_histo->Insert(dabs(Majoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

Primitive_Observable_Base * PT_In_Thrust::Copy() const 
{
  return new PT_In_Thrust(this);
}

//================================================================================
//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(PT_Out_Thrust,PT_Out_Thrust_Getter,"PTOut")

PT_Out_Thrust::PT_Out_Thrust(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{ 
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

PT_Out_Thrust::PT_Out_Thrust(const PT_Out_Thrust * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void PT_Out_Thrust::Evaluate(const Particle_List & pl,double weight,int ncount) 
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    Vec3D Minoraxis = data->Get<Event_Shape_EE_Data>().minoraxis;
    for (Particle_List::const_iterator pit=pl.begin();pit!=pl.end();++pit) {
      p_histo->Insert(dabs(Minoraxis*Vec3D((*pit)->Momentum())),weight,ncount);
    }
  }
}

Primitive_Observable_Base * PT_Out_Thrust::Copy() const 
{
  return new PT_Out_Thrust(this);
}

//================================================================================
//================================================================================
//================================================================================

DEFINE_OBSERVABLE_GETTER(Eta_Thrust,Eta_Thrust_Getter,"EtaThrust")

Eta_Thrust::Eta_Thrust(const String_Matrix & parameters) :
  Event_Shapes_Observable_Base(parameters)
{ 
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable Multi in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"EEShapes";
    m_name="";
    if (m_listname!="EEShapes") m_name=m_listname+string("_");
    m_name+="Multi.dat";
  }
  else {
    if (m_listname=="") m_listname="EEShapes";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="EEShapes") m_name=m_listname+string("_");
      m_name+="Multi";
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "Multiplicity";
      p_histo->SetTitle(title);
    }
  }
}

Eta_Thrust::Eta_Thrust(const Eta_Thrust * old) :
  Event_Shapes_Observable_Base(*old)
{
  m_listname = old->m_listname; // ?
  m_name     = old->m_name;     // ?
}

void Eta_Thrust::Evaluate(const ATOOLS::Blob_List & ,double weight, int ncount)
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
  return new Eta_Thrust(this);
}

