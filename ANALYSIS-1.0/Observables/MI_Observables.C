#include "MI_Observables.H"

using namespace ANALYSIS;

#include "MyStrStream.H"
#include "Kt_Algorithm.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"Charged";
    std::string jetlist=parameters[0].size()>5?parameters[0][5]:"AnalysedJets";
    return new Class(10*(int)(parameters[0][3]=="Log"),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),jetlist,list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list="Charged", jetlist="AnalysedJets", scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="JETLIST") jetlist=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class((scale=="Log")*10,min,max,bins,jetlist,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *const					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|Log [jetlist] [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME)

template <class Class>
Primitive_Observable_Base *const GetOffsetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>5?parameters[0][5]:"Charged";
    std::string jetlist=parameters[0].size()>6?parameters[0][6]:"AnalysedJets";
    return new Class(10*(int)(parameters[0][3]=="Log"),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][4]),jetlist,list);
  }
  else if (parameters.size()<5) return NULL;
  double min=0.0, max=1.0, offset=0.0;
  size_t bins=100;
  std::string list="Charged", jetlist="AnalysedJets", scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="OFFSET") offset=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="JETLIST") jetlist=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class((scale=="Log")*10,min,max,bins,offset,jetlist,list);
}									

#define DEFINE_OFFSET_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *const					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetOffsetObservable<CLASS>(parameters); }

#define DEFINE_OFFSET_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins offset Lin|Log [jetlist] [list]"; }

#define DEFINE_OFFSET_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_OFFSET_GETTER_METHOD(CLASS,NAME);				\
  DEFINE_OFFSET_PRINT_METHOD(NAME)

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
  m_name="MI_Statistics.dat";
  m_type=type;
  m_listname=listname;
  m_splitt_flag=false;
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
				 const int nbins,const std::string &jetlist,
				 const std::string &listname):
  Primitive_Observable_Base(0,detamin,detamax,nbins,NULL)
{
  m_name="FwBwEtaCorr.dat";
  m_listname=listname;
  m_etafw.Initialize(detamin,detamax,nbins);
  m_etafwsq.Initialize(detamin,detamax,nbins);
  m_etafwbw.Initialize(detamin,detamax,nbins);
}

void Forward_Backward_Eta_Correlation::
Evaluate(const ATOOLS::Particle_List &particlelist,double weight,int ncount)
{
  ATOOLS::Histogram etafw(0,m_xmin,m_xmax,m_nbins);
  ATOOLS::Histogram etabw(0,m_xmin,m_xmax,m_nbins);
  for (Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double eta=(*pit)->Momentum().Eta();
    if (eta>0.) etafw.Insert(eta,1.);
    else etabw.Insert(-eta,1.);
  }
  double width=(p_histo->Xmax()-p_histo->Xmin())/p_histo->Nbin();
  for (int i=1;i<=p_histo->Nbin();++i) {
    double eta=p_histo->Xmin()+(i-1)*width;
    double nfw=etafw.Value(i), nbw=etabw.Value(i);
    m_etafw.Add(eta,nfw,ncount);
    m_etafwsq.Add(eta,nfw*nfw,ncount);
    m_etafwbw.Add(eta,nfw*nbw,ncount);
  }
}

Primitive_Observable_Base *Forward_Backward_Eta_Correlation::Copy() const
{
  return new Forward_Backward_Eta_Correlation(m_type,m_xmin,m_xmax,m_nbins,
					      "",m_listname);
}

void Forward_Backward_Eta_Correlation::EndEvaluation(double scale) 
{
  for (size_t i=0;i<(size_t)m_nbins+2;++i) {
    double nfwm=m_etafw.BinContent(i)/m_etafw.BinEntries(i);
    p_histo->Bin((int)i)[0]=
      (m_etafwbw.BinContent(i)/m_etafwbw.BinEntries(i)-nfwm*nfwm)/
      (m_etafwsq.BinContent(i)/m_etafwsq.BinEntries(i)-nfwm*nfwm);
  }
  p_histo->Output();
}

#ifdef SORT
#define SORT_LIST(LISTNAME,PREDICATE)					\
  std::sort(LISTNAME->begin(),LISTNAME->end(),ATOOLS::PREDICATE())	
#else
#define SORT_LIST(LISTNAME,PREDICATE)
#endif

#define TRANSFER_DATA(SCALE)						\
  for (size_t i=0;i<(size_t)m_nbins+2;++i)				\
    m_histogram.Add(histo.BinXMean(i),histo.BinEntries(i)!=0?		\
		    histo.BinContent(i)/histo.BinEntries(i)*SCALE:0)


#define FINISH_HISTOGRAM					\
  for (size_t i=0;i<(size_t)m_nbins+2;++i) {		\
    p_histo->Bin((int)i)[0]=m_histogram.BinEntries(i)!=0?	\
      m_histogram.BinContent(i)/m_histogram.BinEntries(i):0;	\
  }								\
  p_histo->Output()

DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_JetPT,
			 Multiplicity_vs_JetPT_Getter,"NvsJetPT");

Multiplicity_vs_JetPT::Multiplicity_vs_JetPT(const int type,
					     const double ptmin,const double ptmax,
					     const int nbins,const std::string &jetlist,
					     const std::string &listname):
  Primitive_Observable_Base(type,ptmin,ptmax,nbins,NULL),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_PT"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Multiplicity_vs_JetPT::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  m_histogram.Add((*jetlist)[0]->Momentum().PPerp(),weight*particlelist.size());
}
    
Primitive_Observable_Base *Multiplicity_vs_JetPT::Copy() const
{
  return new Multiplicity_vs_JetPT(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

void Multiplicity_vs_JetPT::EndEvaluation(double scale)
{
  FINISH_HISTOGRAM;
}

DEFINE_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_JetPT,
			 Scalar_PT_Sum_vs_JetPT_Getter,"ScPTvsJetPT");

Scalar_PT_Sum_vs_JetPT::
Scalar_PT_Sum_vs_JetPT(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname):
  Primitive_Observable_Base(type,ptmin,ptmax,nbins,NULL),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_PT"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Scalar_PT_Sum_vs_JetPT::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  double pt=0.0;
  for (size_t i=0;i<particlelist.size();++i) pt+=particlelist[i]->Momentum().PPerp();
  m_histogram.Add((*jetlist)[0]->Momentum().PPerp(),weight*pt);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_JetPT::Copy() const
{
  return new Scalar_PT_Sum_vs_JetPT(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

void Scalar_PT_Sum_vs_JetPT::EndEvaluation(double scale)
{
  FINISH_HISTOGRAM;
}

DEFINE_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_JetET,
			 Scalar_PT_Sum_vs_JetET_Getter,"ScPTvsJetET");

Scalar_PT_Sum_vs_JetET::
Scalar_PT_Sum_vs_JetET(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname):
  Primitive_Observable_Base(type,ptmin,ptmax,nbins,NULL),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_ET"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Scalar_PT_Sum_vs_JetET::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_ET);
  double pt=0.0;
  for (size_t i=0;i<particlelist.size();++i) pt+=particlelist[i]->Momentum().PPerp();
  m_histogram.Add((*jetlist)[0]->Momentum().EPerp(),weight*pt);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_JetET::Copy() const
{
  return new Scalar_PT_Sum_vs_JetET(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

void Scalar_PT_Sum_vs_JetET::EndEvaluation(double scale)
{
  FINISH_HISTOGRAM;
}

DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_PT,
			 Multiplicity_vs_PT_Getter,"NvsPT");

Multiplicity_vs_PT::Multiplicity_vs_PT(const int type,
				       const double ptmin,const double ptmax,
				       const int nbins,const std::string &jetlist,
				       const std::string &listname):
  Primitive_Observable_Base(type,ptmin,ptmax,nbins,NULL),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_PT_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Multiplicity_vs_PT::Evaluate(const ATOOLS::Particle_List &particlelist,
				  double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    histo.Add((*pit)->Momentum().PPerp(),weight);
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Multiplicity_vs_PT::Copy() const
{
  return new Multiplicity_vs_PT(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

void Multiplicity_vs_PT::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));			
  FINISH_HISTOGRAM;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Multiplicity_vs_DPhi,
				Multiplicity_vs_DPhi_Getter,"NvsDPhi");

Multiplicity_vs_DPhi::Multiplicity_vs_DPhi(const int type,
					   const double dphimin,const double dphimax,
					   const int nbins,const double offset,
					   const std::string &jetlist,
					   const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist),
  m_offset(offset)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_DPhi_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
  while (m_offset>=360.0) m_offset-=360.0;
  while (m_offset<0.0) m_offset+=360.0;
}
    
void Multiplicity_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double phi=(*pit)->Momentum().DPhi(leadingjet)/M_PI*180.0;
    if (!(phi>0.0) && !(phi<=0.0)) continue;
    histo.Add(phi+m_offset-((int)(phi+m_offset)/360)*360.0,weight);
    phi*=-1.0;
    histo.Add(phi+m_offset-((int)(phi+m_offset)/360)*360.0,weight);
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Multiplicity_vs_DPhi::Copy() const
{
  return new Multiplicity_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlist,m_listname);
}

void Multiplicity_vs_DPhi::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));
  FINISH_HISTOGRAM;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_DPhi,
				Scalar_PT_Sum_vs_DPhi_Getter,"ScPTvsDPhi");

Scalar_PT_Sum_vs_DPhi::Scalar_PT_Sum_vs_DPhi(const int type,
					     const double dphimin,const double dphimax,
					     const int nbins,const double offset,
					     const std::string &jetlist,
					     const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_DPhi_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Scalar_PT_Sum_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double pt=(*pit)->Momentum().PPerp();
    double phi=(*pit)->Momentum().DPhi(leadingjet)/M_PI*180.0;
    if (!(phi>0.0) && !(phi<=0.0)) continue;
    histo.Add(phi+m_offset-((int)(phi+m_offset)/360)*360.0,weight*pt);
    phi*=-1.0;
    histo.Add(phi+m_offset-((int)(phi+m_offset)/360)*360.0,weight*pt);
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_DPhi::Copy() const
{
  return new Scalar_PT_Sum_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				   m_jetlist,m_listname);
}

void Scalar_PT_Sum_vs_DPhi::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));
  FINISH_HISTOGRAM;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_ET_Sum_vs_DPhi,
				Scalar_ET_Sum_vs_DPhi_Getter,"ScETvsDPhi");

Scalar_ET_Sum_vs_DPhi::Scalar_ET_Sum_vs_DPhi(const int type,
					     const double dphimin,const double dphimax,
					     const int nbins,const double offset,
					     const std::string &jetlist,
					     const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SET"+m_listname+"_vs_DPhi_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Scalar_ET_Sum_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double pt=(*pit)->Momentum().PPerp();
    double phi=(*pit)->Momentum().DPhi(leadingjet)/M_PI*180.0;
    if (!(phi>0.0) && !(phi<=0.0)) continue;
    histo.Add(phi+m_offset-((int)(phi+m_offset)/360)*360.0,weight*pt);
    phi*=-1.0;
    histo.Add(phi+m_offset-((int)(phi+m_offset)/360)*360.0,weight*pt);
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Scalar_ET_Sum_vs_DPhi::Copy() const
{
  return new Scalar_ET_Sum_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				   m_jetlist,m_listname);
}

void Scalar_ET_Sum_vs_DPhi::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));
  FINISH_HISTOGRAM;
}


DEFINE_OFFSET_OBSERVABLE_GETTER(Multiplicity_vs_DEta,
				Multiplicity_vs_DEta_Getter,"NvsDEta");

Multiplicity_vs_DEta::Multiplicity_vs_DEta(const int type,
					   const double dphimin,const double dphimax,
					   const int nbins,const double offset,
					   const std::string &jetlist,
					   const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist),
  m_offset(offset)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_DEta_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Multiplicity_vs_DEta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    histo.Add((*pit)->Momentum().DEta(leadingjet)-m_offset,weight);
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Multiplicity_vs_DEta::Copy() const
{
  return new Multiplicity_vs_DEta(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlist,m_listname);
}

void Multiplicity_vs_DEta::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));
  FINISH_HISTOGRAM;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_DEta,
				Scalar_PT_Sum_vs_DEta_Getter,"ScPTvsDEta");

Scalar_PT_Sum_vs_DEta::Scalar_PT_Sum_vs_DEta(const int type,
					   const double dphimin,const double dphimax,
					   const int nbins,const double offset,
					   const std::string &jetlist,
					   const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_DEta_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Scalar_PT_Sum_vs_DEta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    histo.Add((*pit)->Momentum().DEta(leadingjet)-m_offset,
	      weight*(*pit)->Momentum().PPerp());
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_DEta::Copy() const
{
  return new Scalar_PT_Sum_vs_DEta(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlist,m_listname);
}

void Scalar_PT_Sum_vs_DEta::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));
  FINISH_HISTOGRAM;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_ET_Sum_vs_DEta,
				Scalar_ET_Sum_vs_DEta_Getter,"ScETvsDEta");

Scalar_ET_Sum_vs_DEta::Scalar_ET_Sum_vs_DEta(const int type,
					   const double dphimin,const double dphimax,
					   const int nbins,const double offset,
					   const std::string &jetlist,
					   const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SET"+m_listname+"_vs_DEta_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Scalar_ET_Sum_vs_DEta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    histo.Add((*pit)->Momentum().DEta(leadingjet)-m_offset,
	      weight*(*pit)->Momentum().EPerp());
  }
  TRANSFER_DATA(1.0);
}
    
Primitive_Observable_Base *Scalar_ET_Sum_vs_DEta::Copy() const
{
  return new Scalar_ET_Sum_vs_DEta(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlist,m_listname);
}

void Scalar_ET_Sum_vs_DEta::EndEvaluation(double scale)
{
  m_histogram.Scale(m_nbins/(m_xmax-m_xmin));
  FINISH_HISTOGRAM;
}

