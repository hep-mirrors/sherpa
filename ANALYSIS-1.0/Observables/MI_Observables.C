#include "MI_Observables.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

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

#include <iomanip>

DECLARE_GETTER(Transverse_Region_Selector_Getter,"TransReg",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *const 
Transverse_Region_Selector_Getter::operator()(const String_Matrix &parameters) const
{									
  if (parameters.size()<5) return NULL;
  double min=60.0, max=120.0;
  std::string inlist="Charged", outlist="TransCharged", jetlist="AnalysedJets";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="JetList") jetlist=parameters[i][1];
    else if (parameters[i][0]=="PhiMin") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="PhiMax") max=ATOOLS::ToType<double>(parameters[i][1]);
  }
  return new Transverse_Region_Selector(min,max,jetlist,inlist,outlist);
}									

void Transverse_Region_Selector_Getter::
PrintInfo(std::ostream &str,const size_t width) const   
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n" 
     <<std::setw(width+7)<<" "<<"OutList list\n" 
     <<std::setw(width+7)<<" "<<"JetList list\n" 
     <<std::setw(width+7)<<" "<<"PhiMin  min\n" 
     <<std::setw(width+7)<<" "<<"PhiMax  max\n" 
     <<std::setw(width+4)<<" "<<"}"; 
}

DECLARE_GETTER(Leading_PT_Selector_Getter,"JetPTSel",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *const 
Leading_PT_Selector_Getter::operator()(const String_Matrix &parameters) const
{									
  if (parameters.size()<5) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", outlist="LeadJets";
  size_t item=0;
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="PTMin") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="PTMax") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
  }
  return new Leading_PT_Selector(item,min,max,inlist,outlist);
}									

void Leading_PT_Selector_Getter::
PrintInfo(std::ostream &str,const size_t width) const   
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n" 
     <<std::setw(width+7)<<" "<<"OutList list\n" 
     <<std::setw(width+7)<<" "<<"PTMin   min\n" 
     <<std::setw(width+7)<<" "<<"PTMax   max\n" 
     <<std::setw(width+7)<<" "<<"Item    item\n" 
     <<std::setw(width+4)<<" "<<"}"; 
}

DECLARE_GETTER(Leading_ET_Selector_Getter,"JetETSel",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *const 
Leading_ET_Selector_Getter::operator()(const String_Matrix &parameters) const
{									
  if (parameters.size()<5) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", outlist="LeadJets";
  size_t item=0;
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="ETMin") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="ETMax") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
  }
  return new Leading_ET_Selector(item,min,max,inlist,outlist);
}									

void Leading_ET_Selector_Getter::
PrintInfo(std::ostream &str,const size_t width) const   
{ 
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n" 
     <<std::setw(width+7)<<" "<<"OutList list\n" 
     <<std::setw(width+7)<<" "<<"ETMin   min\n" 
     <<std::setw(width+7)<<" "<<"ETMax   max\n" 
     <<std::setw(width+7)<<" "<<"Item    item\n" 
     <<std::setw(width+4)<<" "<<"}"; 
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
  m_splitt_flag=false;
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
  for (size_t i=1;i<=(size_t)p_histo->Nbin();++i) {
    double nfwm=m_etafw.BinContent(i)/m_etafw.BinEntries(i);
    p_histo->Bin((int)i)[0]=
      (m_etafwbw.BinContent(i)/m_etafwbw.BinEntries(i)-nfwm*nfwm)/
      (m_etafwsq.BinContent(i)/m_etafwsq.BinEntries(i)-nfwm*nfwm);
  }
  p_histo->Output();
}

Leading_PT_Selector::
Leading_PT_Selector(const size_t item,const double ptmin,const double ptmax,
		    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(ptmin)+"<Leading_PT<"+
	    ATOOLS::ToString(ptmax)+inlist),
  m_item(item)
{
  m_xmin=ptmin;
  m_xmax=ptmax;
  m_listname=inlist;
}

void Leading_PT_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  if (particlelist.size()<=m_item) return;
  double pt=particlelist[m_item]->Momentum().PPerp();
  if (pt<m_xmin || pt>m_xmax) return;
  outlist->resize(particlelist.size());
  for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
}

Primitive_Observable_Base *Leading_PT_Selector::Copy() const
{
  return new Leading_PT_Selector(m_item,m_xmin,m_xmax,m_listname,m_outlist);
}

void Leading_PT_Selector::EndEvaluation(double scale)
{
}

Leading_ET_Selector::
Leading_ET_Selector(const size_t item,const double etmin,const double etmax,
		    const std::string &inlist,const std::string &outlist):
  m_outlist(outlist!=""?outlist:ATOOLS::ToString(etmin)+"<Leading_ET<"+
	    ATOOLS::ToString(etmax)+inlist),
  m_item(item)
{
  m_xmin=etmin;
  m_xmax=etmax;
  m_listname=inlist;
}

void Leading_ET_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  if (particlelist.size()<=m_item) return;
  double et=particlelist[m_item]->Momentum().EPerp();
  if (et<m_xmin || et>m_xmax) return;
  outlist->resize(particlelist.size());
  for (size_t i=0;i<particlelist.size();++i) 
    (*outlist)[i] = new ATOOLS::Particle(*particlelist[i]);
}

Primitive_Observable_Base *Leading_ET_Selector::Copy() const
{
  return new Leading_ET_Selector(m_item,m_xmin,m_xmax,m_listname,m_outlist);
}

void Leading_ET_Selector::EndEvaluation(double scale)
{
}

Transverse_Region_Selector::
Transverse_Region_Selector(const double phimin,const double phimax,
			   const std::string &jetlist,const std::string &inlist,
			   const std::string &outlist):
  m_jetlist(jetlist),
  m_outlist(outlist!=""?outlist:"Trans"+inlist)
{
  m_xmin=phimin;
  m_xmax=phimax;
  m_listname=inlist;
}

void Transverse_Region_Selector::Evaluate(const ATOOLS::Particle_List &particlelist,
					  double weight,int ncount)
{
  ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
  p_ana->AddParticleList(m_outlist,outlist);
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  m_leadingjet=(*jetlist)[0]->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double phi=(*pit)->Momentum().DPhi(m_leadingjet)/M_PI*180.;
    if (phi>m_xmin && phi<m_xmax) outlist->push_back(new ATOOLS::Particle(**pit));
  }
}

Primitive_Observable_Base *Transverse_Region_Selector::Copy() const
{
  return new Transverse_Region_Selector(m_xmin,m_xmax,m_jetlist,
					m_listname,m_outlist);
}

void Transverse_Region_Selector::EndEvaluation(double scale)
{
}

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
    
DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_JetPT,
			 Multiplicity_vs_JetPT_Getter,"NvsJetPT");

void Multiplicity_vs_JetPT::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  m_histogram.Add((*jetlist)[0]->Momentum().PPerp(),weight*particlelist.size(),ncount);
}
    
Primitive_Observable_Base *Multiplicity_vs_JetPT::Copy() const
{
  return new Multiplicity_vs_JetPT(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

void Multiplicity_vs_JetPT::EndEvaluation(double scale)
{
  // m_histogram.Scale((m_histogram.XMax()-m_histogram.XMin())/m_histogram.Norm());
  for (size_t i=1;i<=(size_t)p_histo->Nbin();++i) {
    p_histo->Bin((int)i)[0]=m_histogram.BinEntries(i)!=0?
      m_histogram.BinContent(i)/m_histogram.BinEntries(i):0;
  }
  p_histo->Output();
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
}
    
void Multiplicity_vs_PT::Evaluate(const ATOOLS::Particle_List &particlelist,
				     double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  // ATOOLS::Vec4D leadingjet=(*jetlist->begin())->Momentum();
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    p_histo->Insert((*pit)->Momentum().PPerp(),weight,ncount);
  }
}
    
Primitive_Observable_Base *Multiplicity_vs_PT::Copy() const
{
  return new Multiplicity_vs_PT(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Multiplicity_vs_DPhi,
			 Multiplicity_vs_DPhi_Getter,"NvsDPhi");

Multiplicity_vs_DPhi::Multiplicity_vs_DPhi(const int type,
					   const double dphimin,const double dphimax,
					   const int nbins,const std::string &jetlist,
					   const std::string &listname):
  Primitive_Observable_Base(type,dphimin,dphimax,nbins,NULL),
  m_jetlist(jetlist)
{
  m_listname=listname;
  m_name="N"+m_listname+"_vs_DPhi_"+m_jetlist+".dat";
  m_histogram.Initialize(m_xmin,m_xmax,m_nbins);
}
    
void Multiplicity_vs_DPhi::Evaluate(const ATOOLS::Particle_List &particlelist,
				    double weight,int ncount)
{
  ATOOLS::Particle_List *jetlist=p_ana->GetParticleList(m_jetlist);
  if (jetlist->size()==0) return;
  ATOOLS::Vec4D leadingjet=(*jetlist)[0]->Momentum();
//   ATOOLS::Vec3D perp=ATOOLS::cross(leadingjet,ATOOLS::Vec3D::ZVEC);
//   if (jetlist->size()>1) {
//     ATOOLS::Vec3D test=ATOOLS::cross(leadingjet,(*jetlist)[1]->Momentum());
//     test=ATOOLS::cross(test,leadingjet);
//     if (test.Abs()>0.0) perp=test;
//   }
//   else if (perp.Abs()==0.0) {
//     perp=ATOOLS::cross(leadingjet,ATOOLS::Vec3D::XVEC);
//   }
  AMISIC::Amisic_Histogram<double> histo;
  histo.Initialize(m_xmin,m_xmax,m_nbins);
  for (ATOOLS::Particle_List::const_iterator pit=particlelist.begin();
       pit!=particlelist.end();++pit) {
    double phi=(*pit)->Momentum().DPhi(leadingjet);
//     if (perp*(ATOOLS::Vec3D)(*pit)->Momentum()<0.0) phi*=-1.0;
//     phi+=1.5*M_PI;
//     while (phi>=2.0*M_PI) phi-=2.0*M_PI;
//     while (phi<0.0) phi+=2.0*M_PI;
    phi-=0.5*M_PI;
    histo.Add((phi>=0.0?phi:2.0*M_PI+phi)/M_PI*180.,weight);
    phi=M_PI-phi;
    histo.Add(phi/M_PI*180.,weight);
  }
  for (size_t i=1;i<=(size_t)p_histo->Nbin()-1;++i) {
    m_histogram.Add(histo.BinXMean(i),histo.BinContent(i),ncount);
  }
}
    
Primitive_Observable_Base *Multiplicity_vs_DPhi::Copy() const
{
  return new Multiplicity_vs_DPhi(m_type,m_xmin,m_xmax,m_nbins,m_jetlist,m_listname);
}

void Multiplicity_vs_DPhi::EndEvaluation(double scale)
{
  for (size_t i=1;i<=(size_t)p_histo->Nbin();++i) {
    p_histo->Bin((int)i)[0]=m_histogram.BinEntries(i)!=0?
      m_histogram.BinContent(i)/m_histogram.BinEntries(i):0;
  }
  p_histo->Output();
}

