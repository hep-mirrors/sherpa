#include "MI_Observables.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class,class Getter>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"Analysed";
    return new Class(10*(int)(parameters[0][4]=="Log"),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100, scale=0;
  std::string list="Analysed";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(scale,min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME,TAG)				\
  Primitive_Observable_Base *const					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS,NAME>(parameters); }

#define DEFINE_PRINT_METHOD(CLASS,NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|Log [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME,TAG);					\
  DEFINE_PRINT_METHOD(CLASS,NAME)

DECLARE_GETTER(MI_Statistics_Getter,"MIStats",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *const 
MI_Statistics_Getter::operator()(const String_Matrix &parameters) const
{
  std::string listname="Analysed";
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new MI_Statistics(listname);
}

void MI_Statistics_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

#include "Primitive_Analysis.H"

#include <fstream>

using namespace ATOOLS;

MI_Statistics::MI_Statistics(const std::string & listname, int type):
  Primitive_Observable_Base(type,0,100,100,NULL) 
{
  m_name  = "MI_Statistics.dat";
  m_type  = type;
  m_listname    = listname;
  m_splitt_flag = false;
}

void MI_Statistics::Evaluate(const Blob_List &  blobs,double weight,int ncount)
{
  unsigned int number=0;
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type()==btp::Hard_Collision) {
      ++number;
    }
  }
  p_histo->Insert(number,weight,ncount);
}

Primitive_Observable_Base * MI_Statistics::Copy() const 
{
  return new MI_Statistics(m_listname,m_type);
}

DEFINE_OBSERVABLE_GETTER(Forward_Backward_Eta_Correlation,
			 Forward_Backward_Eta_Correlation_Getter,"EtaCorr");

Forward_Backward_Eta_Correlation::
Forward_Backward_Eta_Correlation(const int type,
				 const double detamin,const double detamax,
				 const int nbins,const std::string &listname):
  Primitive_Observable_Base(type,detamin,detamax,nbins,NULL),
  m_etafw(0,detamin/2.,detamax/2.+1.,nbins+1),
  m_etabw(0,detamin/2.,detamax/2.+1.,nbins+1)
{
  m_name="Eta_Correlator.dat";
  m_type=type;
  m_listname=listname;
  m_splitt_flag=false;
}

void Forward_Backward_Eta_Correlation::
Evaluate(const ATOOLS::Particle_List &particles,double weight,int ncount)
{
  for (Particle_List::const_iterator pit=particles.begin();
       pit!=particles.end();++pit) {
    if (!(*pit)->Flav().Charge()!=0.) continue;
    double eta=(*pit)->Momentum().Eta();
    if (eta>0.) m_etafw.Insert(eta,1.);
    else m_etabw.Insert(eta,1.);
  }
  m_etafw.Reset();
  m_etabw.Reset();
}

Primitive_Observable_Base *Forward_Backward_Eta_Correlation::Copy() const
{
  return new Forward_Backward_Eta_Correlation(m_type,m_xmin,m_xmax,m_nbins,NULL);
}

