#include "AddOns/Analysis/Observables/MI_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Triggers/Kt_Algorithm.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <algorithm>

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:"Charged";
    std::string jetlist=parameters[0].size()>5?parameters[0][5]:"AnalysedJets";
    return new Class(HistogramType(parameters[0][3]),
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
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Bins") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Scale") scale=parameters[i][1];
    else if (parameters[i][0]=="JetList") jetlist=parameters[i][1];
    else if (parameters[i][0]=="List") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,jetlist,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					        \
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [jetlist] [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

template <class Class>
Primitive_Observable_Base *const GetOffsetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>5?parameters[0][5]:"Charged";
    std::string jetlist=parameters[0].size()>6?parameters[0][6]:"AnalysedJets";
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     std::vector<std::string>(1,jetlist),list);
  }
  else if (parameters.size()<5) return NULL;
  double min=0.0, max=1.0, offset=0.0;
  size_t bins=100;
  std::string list="Charged", scale="Lin";
  std::vector<std::string> jetlists;
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Bins") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Offset") offset=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Scale") scale=parameters[i][1];
    else if (parameters[i][0]=="JetList") jetlists.push_back(parameters[i][1]);
    else if (parameters[i][0]=="List") list=parameters[i][1];
  }
  if (jetlists.size()==0) jetlists.push_back("Jets");
  return new Class(HistogramType(scale),min,max,bins,offset,jetlists,list);
}									

#define DEFINE_OFFSET_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					        \
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetOffsetObservable<CLASS>(parameters); }

#define DEFINE_OFFSET_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins offset Lin|LinErr|Log|LogErr [jetlist] [list]"; }

#define DEFINE_OFFSET_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_OFFSET_GETTER_METHOD(CLASS,NAME)				\
  DEFINE_OFFSET_PRINT_METHOD(NAME)

DECLARE_GETTER(MI_Statistics_Getter,"MIStats",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base * 
MI_Statistics_Getter::operator()(const Argument_Matrix &parameters) const
{
  size_t scales=5;
  std::string listname=finalstate_list;
  if (parameters.size()>0 && parameters[0].size()>0) {
    scales=ATOOLS::ToType<int>(parameters[0][0]);
    if (parameters[0].size()>1) listname=parameters[0][1];
  }
  return new MI_Statistics(ATOOLS::Max((size_t)1,scales),listname);
}

void MI_Statistics_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[scales] [list]"; 
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <fstream>

using namespace ATOOLS;

MI_Statistics::MI_Statistics(const size_t scales,const std::string & listname, 
			     int type):
  Primitive_Observable_Base(type,0,100,100) 
{
  m_name="MI_Statistics.dat";
  m_type=type;
  m_listname=listname;
  m_scales.resize(scales);
  double max=ATOOLS::rpa.gen.Ecms()/2.0;
  for (size_t i=0;i<scales;++i) 
    m_scales[i] = new ATOOLS::Histogram(10,max*1.0e-5,max,200);
}

MI_Statistics::~MI_Statistics()
{
  for (size_t i=0;i<m_scales.size();++i) delete m_scales[i];
}

void MI_Statistics::Evaluate(const Blob_List &  blobs,double weight,double ncount)
{
  double hard=0.0;
  unsigned int number=0;
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    double scale=0.0;
    if ((*bit)->Type()==ATOOLS::btp::Hard_Collision) {
      ++number;
      ATOOLS::Blob_Data_Base *info=(*(*bit))["MI_Scale"];
      if (info!=NULL) scale=info->Get<double>();
    }
    else if ((*bit)->Type()==ATOOLS::btp::Signal_Process) {
      ATOOLS::Blob_Data_Base *info=(*(*bit))["MI_Scale"];
      if (info!=NULL) hard=info->Get<double>();
    }
    if (number-1<m_scales.size() && scale!=0.0) {
      m_scales[number-1]->Insert(scale,weight,ncount);
    }
  }
  if (number==0) m_scales[0]->Insert(hard,weight,ncount);
  p_histo->Insert((double)number,weight,ncount);
}

Primitive_Observable_Base &MI_Statistics::
operator+=(const Primitive_Observable_Base &obs)
{
  const MI_Statistics *mi=dynamic_cast<const MI_Statistics*>(&obs);
  if (mi!=NULL) {
    (*p_histo)+=*mi->p_histo;
    for (size_t i=0;i<m_scales.size();++i) (*m_scales[i])+=*mi->m_scales[i];
  }
  return *this;
}

void MI_Statistics::EndEvaluation(double scale) 
{
  p_histo->Finalize();
  if (scale!=1.) p_histo->Scale(scale);
  for (size_t i=0;i<m_scales.size();++i) {
    m_scales[i]->Finalize();
    if (scale!=1.) m_scales[i]->Scale(scale);
  }
}

void MI_Statistics::Restore(double scale) 
{
  if (scale!=1.) p_histo->Scale(scale);
  p_histo->Restore();
  for (size_t i=0;i<m_scales.size();++i) {
    if (scale!=1.) m_scales[i]->Scale(scale);
    m_scales[i]->Restore();
  }
}

void MI_Statistics::Output(const std::string & pname) 
{
  if (p_histo) {
    ATOOLS::MakeDir(pname); 
    p_histo->Output((pname+std::string("/")+m_name).c_str());
  }
  for (size_t i=0;i<m_scales.size();++i) {
    m_scales[i]->Output(pname+"/MIScale"+m_listname+"_"+ATOOLS::ToString(i)+".dat");
  }
}

Primitive_Observable_Base * MI_Statistics::Copy() const 
{
  return new MI_Statistics(m_scales.size(),m_listname,m_type);
}

#define SORT
#ifdef SORT
#define SORT_LIST(LISTNAME,PREDICATE)					\
  std::sort(LISTNAME->begin(),LISTNAME->end(),PREDICATE())	
#else
#define SORT_LIST(LISTNAME,PREDICATE)
#endif

DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_JetPT,
			 Multiplicity_vs_JetPT_Getter,"NvsJetPT")

Multiplicity_vs_JetPT::
Multiplicity_vs_JetPT(const int type,
		      const double ptmin,const double ptmax,
		      const int nbins,const std::string &jetlist,
		      const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_PT"+m_jetlist+".dat";
}
    
void Multiplicity_vs_JetPT::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  double ptjet=(*jetlist)[0]->Momentum().PPerp();
  p_obs->Insert(ptjet,particlelist.size()*weight);
  p_norm->Insert(ptjet,weight);
}
    
Primitive_Observable_Base *Multiplicity_vs_JetPT::Copy() const
{
  Multiplicity_vs_JetPT *obs = 
    new Multiplicity_vs_JetPT(m_type,m_xmin,m_xmax,m_nbins,
			      m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}


DEFINE_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_JetPT,
			 Scalar_PT_Sum_vs_JetPT_Getter,"ScPTvsJetPT")

Scalar_PT_Sum_vs_JetPT::
Scalar_PT_Sum_vs_JetPT(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_PT"+m_jetlist+".dat";
}
    
void Scalar_PT_Sum_vs_JetPT::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  double pt=0.0;
  for (size_t i=0;i<particlelist.size();++i) 
    pt+=particlelist[i]->Momentum().PPerp();
  double ptjet=(*jetlist)[0]->Momentum().PPerp();
  p_obs->Insert(ptjet,pt*weight);
  p_norm->Insert(ptjet,weight);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_JetPT::Copy() const
{
  Scalar_PT_Sum_vs_JetPT *obs = 
    new Scalar_PT_Sum_vs_JetPT(m_type,m_xmin,m_xmax,m_nbins,
			      m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_JetET,
			 Scalar_PT_Sum_vs_JetET_Getter,"ScPTvsJetET")

Scalar_PT_Sum_vs_JetET::
Scalar_PT_Sum_vs_JetET(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_ET"+m_jetlist+".dat";
}
    
void Scalar_PT_Sum_vs_JetET::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_ET);
  double pt=0.0;
  for (size_t i=0;i<particlelist.size();++i) 
    pt+=particlelist[i]->Momentum().PPerp();
  double etjet=(*jetlist)[0]->Momentum().EPerp();
  p_obs->Insert(etjet,pt*weight);
  p_norm->Insert(etjet,weight);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_JetET::Copy() const
{
  Scalar_PT_Sum_vs_JetET *obs = 
    new Scalar_PT_Sum_vs_JetET(m_type,m_xmin,m_xmax,m_nbins,
			      m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_PT,
			 Multiplicity_vs_PT_Getter,"NvsPT")

Multiplicity_vs_PT::Multiplicity_vs_PT(const int type,
				       const double ptmin,const double ptmax,
				       const int nbins,
				       const std::string &jetlist,
				       const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_PT_"+m_jetlist+".dat";
}
    
void Multiplicity_vs_PT::Evaluate(const ATOOLS::Particle_List &particlelist,
				  double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double pt=(*pit)->Momentum().PPerp();
    p_obs->Insert(pt,weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Multiplicity_vs_PT::Copy() const
{
  Multiplicity_vs_PT *obs = 
    new Multiplicity_vs_PT(m_type,m_xmin,m_xmax,m_nbins,
			   m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Multiplicity_vs_DPhi,
				Multiplicity_vs_DPhi_Getter,"NvsDPhi")

Multiplicity_vs_DPhi::
Multiplicity_vs_DPhi(const int type,
		     const double dphimin,const double dphimax,
		     const int nbins,const double offset,
		     const std::vector<std::string> &jetlists,
		     const std::string &listname):
  Normalized_Observable(type,dphimin,dphimax,nbins),
  m_jetlists(jetlists.empty()?std::vector<std::string>(1,"Jets"):jetlists),
  m_offset(offset)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_DPhi_"+m_jetlists[0]+".dat";
  while (m_offset>=360.0) m_offset-=360.0;
  while (m_offset<0.0) m_offset+=360.0;
}
    
void Multiplicity_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist1=p_ana->GetParticleList(m_jetlists[0]);
  if (jetlist1->empty()) return;
  ATOOLS::Particle_List *jetlist2=NULL;
  ATOOLS::Vec4D leadingjet2;
  if (m_jetlists.size()>1) {
    if ((jetlist2=p_ana->GetParticleList(m_jetlists[1]))->empty()) return;
    SORT_LIST(jetlist2,Order_PT);
    leadingjet2=(*jetlist2)[0]->Momentum();
  }
  SORT_LIST(jetlist1,Order_PT);
  ATOOLS::Vec4D leadingjet1=(*jetlist1)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double dphi1=(*pit)->Momentum().DPhi(leadingjet1)/M_PI*180.0;
    if (!(dphi1>0.0) && !(dphi1<=0.0)) continue;
    double cosdphi2=jetlist2==NULL?0.0:(*pit)->Momentum().CosDPhi(leadingjet2);
    if (cosdphi2>=0.0) {
      double cur=dphi1+m_offset-((int)(dphi1+m_offset)/360)*360.0;
      p_obs->Insert(cur,weight);
    }
    if (cosdphi2<=0.0) {
      dphi1*=-1.0;
      double cur=dphi1+m_offset-((int)(dphi1+m_offset)/360)*360.0;
      p_obs->Insert(cur,weight);
    }
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Multiplicity_vs_DPhi::Copy() const
{
  Multiplicity_vs_DPhi *obs = 
    new Multiplicity_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlists,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_DPhi,
				Scalar_PT_Sum_vs_DPhi_Getter,"ScPTvsDPhi")

Scalar_PT_Sum_vs_DPhi::
Scalar_PT_Sum_vs_DPhi(const int type,
		      const double dphimin,const double dphimax,
		      const int nbins,const double offset,
		      const std::vector<std::string> &jetlists,
		      const std::string &listname):
  Normalized_Observable(type,dphimin,dphimax,nbins),
  m_jetlists(jetlists.empty()?std::vector<std::string>(1,"Jets"):jetlists),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_DPhi_"+m_jetlists[0]+".dat";
}
    
void Scalar_PT_Sum_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist1=p_ana->GetParticleList(m_jetlists[0]);
  if (jetlist1->empty()) return;
  ATOOLS::Particle_List *jetlist2=NULL;
  ATOOLS::Vec4D leadingjet2;
  if (m_jetlists.size()>1) {
    if ((jetlist2=p_ana->GetParticleList(m_jetlists[1]))->empty()) return;
    SORT_LIST(jetlist2,Order_PT);
    leadingjet2=(*jetlist2)[0]->Momentum();
  }
  SORT_LIST(jetlist1,Order_PT);
  ATOOLS::Vec4D leadingjet1=(*jetlist1)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double pt=(*pit)->Momentum().PPerp();
    double dphi1=(*pit)->Momentum().DPhi(leadingjet1)/M_PI*180.0;
    if (!(dphi1>0.0) && !(dphi1<=0.0)) continue;
    double cosdphi2=jetlist2==NULL?0.0:(*pit)->Momentum().CosDPhi(leadingjet2);
    if (cosdphi2>=0.0) {
      double cur=dphi1+m_offset-((int)(dphi1+m_offset)/360)*360.0;
      p_obs->Insert(cur,pt*weight);
    }
    if (cosdphi2<=0.0) {
      dphi1*=-1.0;
      double cur=dphi1+m_offset-((int)(dphi1+m_offset)/360)*360.0;
      p_obs->Insert(cur,pt*weight);
    }
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_DPhi::Copy() const
{
  Scalar_PT_Sum_vs_DPhi *obs = 
    new Scalar_PT_Sum_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlists,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_ET_Sum_vs_DPhi,
				Scalar_ET_Sum_vs_DPhi_Getter,"ScETvsDPhi")

Scalar_ET_Sum_vs_DPhi::
Scalar_ET_Sum_vs_DPhi(const int type,
		      const double dphimin,const double dphimax,
		      const int nbins,const double offset,
		      const std::vector<std::string> &jetlists,
		      const std::string &listname):
  Normalized_Observable(type,dphimin,dphimax,nbins),
  m_jetlists(jetlists.empty()?std::vector<std::string>(1,"Jets"):jetlists),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SET"+m_listname+"_vs_DPhi_"+m_jetlists[0]+".dat";
}
    
void Scalar_ET_Sum_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist1=p_ana->GetParticleList(m_jetlists[0]);
  if (jetlist1->empty()) return;
  ATOOLS::Particle_List *jetlist2=NULL;
  ATOOLS::Vec4D leadingjet2;
  if (m_jetlists.size()>1) {
    if ((jetlist2=p_ana->GetParticleList(m_jetlists[1]))->empty()) return;
    SORT_LIST(jetlist2,Order_PT);
    leadingjet2=(*jetlist2)[0]->Momentum();
  }
  SORT_LIST(jetlist1,Order_PT);
  ATOOLS::Vec4D leadingjet1=(*jetlist1)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double et=(*pit)->Momentum().EPerp();
    double dphi1=(*pit)->Momentum().DPhi(leadingjet1)/M_PI*180.0;
    if (!(dphi1>0.0) && !(dphi1<=0.0)) continue;
    double cosdphi2=jetlist2==NULL?0.0:(*pit)->Momentum().CosDPhi(leadingjet2);
    if (cosdphi2>=0.0) {
      double cur=dphi1+m_offset-((int)(dphi1+m_offset)/360)*360.0;
      p_obs->Insert(cur,et*weight);
    }
    if (cosdphi2<=0.0) {
      dphi1*=-1.0;
      double cur=dphi1+m_offset-((int)(dphi1+m_offset)/360)*360.0;
      p_obs->Insert(cur,et*weight);
    }
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Scalar_ET_Sum_vs_DPhi::Copy() const
{
  Scalar_ET_Sum_vs_DPhi *obs = 
    new Scalar_ET_Sum_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_offset,
				  m_jetlists,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Multiplicity_vs_DEta,
				Multiplicity_vs_DEta_Getter,"NvsDEta")

Multiplicity_vs_DEta::
Multiplicity_vs_DEta(const int type,
		     const double dphimin,const double dphimax,
		     const int nbins,const double offset,
		     const std::vector<std::string> &jetlists,
		     const std::string &listname):
  Normalized_Observable(type,dphimin,dphimax,nbins),
  m_jetlist(jetlists.empty()?"":jetlists[0]),
  m_offset(offset)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_DEta_"+m_jetlist+".dat";
}
    
void Multiplicity_vs_DEta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double cur=(*pit)->Momentum().DEta(leadingjet)-m_offset;
    p_obs->Insert(cur,weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Multiplicity_vs_DEta::Copy() const
{
  Multiplicity_vs_DEta *obs = 
    new Multiplicity_vs_DEta(m_type,m_xmin,m_xmax,m_nbins,m_offset,
			     std::vector<std::string>(1,m_jetlist),
			     m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_DEta,
				Scalar_PT_Sum_vs_DEta_Getter,"ScPTvsDEta")

Scalar_PT_Sum_vs_DEta::
Scalar_PT_Sum_vs_DEta(const int type,
		      const double dphimin,const double dphimax,
		      const int nbins,const double offset,
		      const std::vector<std::string> &jetlists,
		      const std::string &listname):
  Normalized_Observable(type,dphimin,dphimax,nbins),
  m_jetlist(jetlists.empty()?"":jetlists[0]),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_DEta_"+m_jetlist+".dat";
}
    
void Scalar_PT_Sum_vs_DEta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double cur=(*pit)->Momentum().DEta(leadingjet)-m_offset;
    p_obs->Insert(cur,(*pit)->Momentum().PPerp()*weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_DEta::Copy() const
{
  Scalar_PT_Sum_vs_DEta *obs = 
    new Scalar_PT_Sum_vs_DEta(m_type,m_xmin,m_xmax,m_nbins,m_offset,
			     std::vector<std::string>(1,m_jetlist),
			     m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OFFSET_OBSERVABLE_GETTER(Scalar_ET_Sum_vs_DEta,
				Scalar_ET_Sum_vs_DEta_Getter,"ScETvsDEta")

Scalar_ET_Sum_vs_DEta::
Scalar_ET_Sum_vs_DEta(const int type,
		      const double dphimin,const double dphimax,
		      const int nbins,const double offset,
		      const std::vector<std::string> &jetlists,
		      const std::string &listname):
  Normalized_Observable(type,dphimin,dphimax,nbins),
  m_jetlist(jetlists.empty()?"":jetlists[0]),
  m_offset(offset)
{
  m_listname=listname;
  m_name="SET"+m_listname+"_vs_DEta_"+m_jetlist+".dat";
}
    
void Scalar_ET_Sum_vs_DEta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double cur=(*pit)->Momentum().DEta(leadingjet)-m_offset;
    p_obs->Insert(cur,(*pit)->Momentum().EPerp()*weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Scalar_ET_Sum_vs_DEta::Copy() const
{
  Scalar_ET_Sum_vs_DEta *obs = 
    new Scalar_ET_Sum_vs_DEta(m_type,m_xmin,m_xmax,m_nbins,m_offset,
			     std::vector<std::string>(1,m_jetlist),
			     m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_Eta,
			 Multiplicity_vs_Eta_Getter,"NvsEta")

Multiplicity_vs_Eta::
Multiplicity_vs_Eta(const int type,
		      const double ptmin,const double ptmax,
		      const int nbins,const std::string &jetlist,
		      const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_Eta_"+m_jetlist+".dat";
}
    
void Multiplicity_vs_Eta::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double cur=(*pit)->Momentum().Eta();
    p_obs->Insert(cur,weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Multiplicity_vs_Eta::Copy() const
{
  Multiplicity_vs_Eta *obs = 
    new Multiplicity_vs_Eta(m_type,m_xmin,m_xmax,m_nbins,
			    m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OBSERVABLE_GETTER(Scalar_PT_Sum_vs_Eta,
			 Scalar_PT_Sum_vs_Eta_Getter,"ScPTvsEta")

Scalar_PT_Sum_vs_Eta::
Scalar_PT_Sum_vs_Eta(const int type,
		     const double ptmin,const double ptmax,
		     const int nbins,const std::string &jetlist,
		     const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="SPT"+m_listname+"_vs_Eta_"+m_jetlist+".dat";
}
    
void Scalar_PT_Sum_vs_Eta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double cur=(*pit)->Momentum().Eta();
    p_obs->Insert(cur,(*pit)->Momentum().PPerp()*weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Scalar_PT_Sum_vs_Eta::Copy() const
{
  Scalar_PT_Sum_vs_Eta *obs = 
    new Scalar_PT_Sum_vs_Eta(m_type,m_xmin,m_xmax,m_nbins,
			     m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OBSERVABLE_GETTER(Scalar_ET_Sum_vs_Eta,
			 Scalar_ET_Sum_vs_Eta_Getter,"ScETvsEta")

Scalar_ET_Sum_vs_Eta::
Scalar_ET_Sum_vs_Eta(const int type,
		     const double ptmin,const double ptmax,
		     const int nbins,const std::string &jetlist,
		     const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="SET"+m_listname+"_vs_Eta_"+m_jetlist+".dat";
}
    
void Scalar_ET_Sum_vs_Eta::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double cur=(*pit)->Momentum().Eta();
    p_obs->Insert(cur,(*pit)->Momentum().EPerp()*weight);
  }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *Scalar_ET_Sum_vs_Eta::Copy() const
{
  Scalar_ET_Sum_vs_Eta *obs = 
    new Scalar_ET_Sum_vs_Eta(m_type,m_xmin,m_xmax,m_nbins,
			     m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}

DEFINE_OBSERVABLE_GETTER(MIScale_vs_JetPT,
			 MIScale_vs_JetPT_Getter,"MISvsJetPT")

MIScale_vs_JetPT::
MIScale_vs_JetPT(const int type,
		      const double ptmin,const double ptmax,
		      const int nbins,const std::string &jetlist,
		      const std::string &listname):
  Normalized_Observable(type,ptmin,ptmax,nbins),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="MIScale"+m_listname+"_vs_PT"+m_jetlist+".dat";
}
    
void MIScale_vs_JetPT::Evaluate(const ATOOLS::Blob_List &bloblist, 
				double weight,double ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  SORT_LIST(jetlist,Order_PT);
  double ptjet=(*jetlist)[0]->Momentum().PPerp();
  for (ATOOLS::Blob_List::const_iterator bit=bloblist.begin();
       bit!=bloblist.end();++bit)
    if ((*bit)->Type()==ATOOLS::btp::Hard_Collision) {
      ATOOLS::Blob_Data_Base *info=(*(*bit))["MI_Scale"];
      if (info!=NULL) {
	double scale=info->Get<double>();
	p_obs->Insert(ptjet,scale*weight);
      }
      break;
    }
  for (int i=0;i<m_nbins+2;++i) p_norm->Insert(i,weight);
}
    
Primitive_Observable_Base *MIScale_vs_JetPT::Copy() const
{
  MIScale_vs_JetPT *obs = 
    new MIScale_vs_JetPT(m_type,m_xmin,m_xmax,m_nbins,
			 m_jetlist,m_listname);
  obs->m_copied=true;
  return obs;
}
