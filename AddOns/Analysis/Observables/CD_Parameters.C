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
	       Analysis_Object,Argument_Matrix);

Analysis_Object * ATOOLS::Getter<Analysis_Object,Argument_Matrix,
				 CD_Parameter_Calculator>::operator()(const Argument_Matrix &parameters) const
{
  std::string listname=finalstate_list;
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new CD_Parameter_Calculator(listname);
}

void ATOOLS::Getter<Analysis_Object,Argument_Matrix,CD_Parameter_Calculator>::PrintInfo(std::ostream &str,const size_t width) const
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
  ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,CLASS>::operator()(const Argument_Matrix &parameters) const \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list] -> CDCalc"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Argument_Matrix);	\
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
