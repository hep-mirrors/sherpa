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
    return new Class(10*(int)(parameters[0][3]=="Log"),
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
    else if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class((scale=="Log")*10,min,max,bins,list);
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

#include "Primitive_Analysis.H"

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
			 Forward_Backward_Eta_Correlation_Getter,"FwBwEta");

Forward_Backward_Eta_Correlation::
Forward_Backward_Eta_Correlation(const int type,
				 const double detamin,const double detamax,
				 const int nbins,const std::string &listname):
  Primitive_Observable_Base(0,detamin,detamax,nbins,NULL)
{
  m_name="FwBwEtaCorr.dat";
  m_listname=listname;
  m_splitt_flag=false;
  m_etafw.Initialize(detamin,detamax,nbins);
  m_etafwsq.Initialize(detamin,detamax,nbins);
  m_etafwbw.Initialize(detamin,detamax,nbins);
}

void Forward_Backward_Eta_Correlation::
Evaluate(const ATOOLS::Blob_List &bloblist,double weight,int ncount)
{
  ATOOLS::Histogram etafw(0,m_xmin,m_xmax,m_nbins);
  ATOOLS::Histogram etabw(0,m_xmin,m_xmax,m_nbins);
  ATOOLS::Particle_List *list=p_ana->GetParticleList(m_listname);
  for (Particle_List::const_iterator pit=list->begin();
       pit!=list->end();++pit) {
    double eta=(*pit)->Momentum().Eta();
    if (eta>0.) etafw.Insert(eta,1.);
    else etabw.Insert(-eta,1.);
  }
  double width=(p_histo->Xmax()-p_histo->Xmin())/p_histo->Nbin();
  for (int i=1;i<=p_histo->Nbin();++i) {
    double eta=p_histo->Xmin()+(i-1)*width;
    double nfw=etafw.Value(i), nbw=etabw.Value(i);
    m_etafw.Add(eta,nfw);
    m_etafwsq.Add(eta,nfw*nfw);
    m_etafwbw.Add(eta,nfw*nbw);
  }
}

Primitive_Observable_Base *Forward_Backward_Eta_Correlation::Copy() const
{
  return new Forward_Backward_Eta_Correlation(m_type,m_xmin,m_xmax,m_nbins,
					      m_listname);
}

void Forward_Backward_Eta_Correlation::EndEvaluation(double scale) 
{
  for (unsigned int i=1;i<=(unsigned int)p_histo->Nbin();++i) {
    double nfwm=m_etafw.BinContent(i)/m_etafw.BinEntries(i);
    p_histo->Bin((int)i)[0]=
      (m_etafwbw.BinContent(i)/m_etafwbw.BinEntries(i)-nfwm*nfwm)/
      (m_etafwsq.BinContent(i)/m_etafwsq.BinEntries(i)-nfwm*nfwm);
  }
  p_histo->Output();
}

