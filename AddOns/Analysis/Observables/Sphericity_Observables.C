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

DECLARE_GETTER(Sphericity_Calculator_Getter,"SphCalc",
	       Analysis_Object,Argument_Matrix);

Analysis_Object * 
Sphericity_Calculator_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string listname=finalstate_list;
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new Sphericity_Calculator(listname);
}

void Sphericity_Calculator_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:finalstate_list;
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list] -> SphCalc"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)


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
  double lambda1=0., lambda2=0., lambda3=0.;
  if (pl->size()>0) {
    m_tensor.Calculate(*pl);
    lambda1=m_tensor.EigenValue(0);
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
