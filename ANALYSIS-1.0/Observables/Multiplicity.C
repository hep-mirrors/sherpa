#include "Multiplicity.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"Analysed";
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list="Analysed", scale="Lin";
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
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME)

#include "MathTools.H"

using namespace ATOOLS;

DEFINE_OBSERVABLE_GETTER(Multiplicity,Multiplicity_Getter,"Multi");

Multiplicity::Multiplicity(int type,double xmin,double xmax,int nbins,
			   const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL)
{
  if (listname!="") {
    m_listname = listname;
    m_name = listname+"_multi.dat";
  }
  else
    m_name = "multi.dat";
}
 
void Multiplicity::Evaluate(const ATOOLS::Particle_List & pl,
			    double weight, int ncount)
{
  p_histo->Insert(pl.size(),weight,ncount); 
}


Primitive_Observable_Base * Multiplicity::Copy() const {
  return new Multiplicity(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}




