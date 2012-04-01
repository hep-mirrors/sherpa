#include "AddOns/Analysis/Observables/Jet_Observables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    std::string list=parameters[0].size()>7?parameters[0][7]:finalstate_list;
    return new Class(HistogramType(parameters[0][6]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),
		     ATOOLS::ToType<int>(parameters[0][5]),list);
  }
  else if (parameters.size()<7) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100, nmin=1, nmax=10, mode=1;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="MODE") mode=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NMIN") nmin=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NMAX") nmax=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,mode,nmin,nmax,list);
}									

template <>
Primitive_Observable_Base *GetObservable<Jet_Differential_Rates>(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    std::string list=parameters[0].size()>7?parameters[0][7]:finalstate_list;
    std::string reflist=parameters[0].size()>8?parameters[0][8]:"";
    return new Jet_Differential_Rates(HistogramType(parameters[0][6]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<int>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),
		     ATOOLS::ToType<int>(parameters[0][5]),list,reflist);
  }
  else if (parameters.size()<7) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100, nmin=1, nmax=10, mode=1;
  std::string list=finalstate_list, scale="Lin";
  std::string reflist="";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="MODE") mode=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NMIN") nmin=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="NMAX") nmax=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
    else if (parameters[i][0]=="REF")  reflist=parameters[i][1];
  }
  return new Jet_Differential_Rates(HistogramType(scale),min,max,bins,mode,nmin,nmax,list,reflist);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins mode nmin nmax Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;

Jet_Observable_Base::Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), m_mode(mode), m_minn(minn), m_maxn(maxn)
{
  m_listname=listname;
  m_name  = std::string("jet_");
  if (listname!=finalstate_list) m_name=listname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }
  if (p_histo) delete p_histo; 
  p_histo =  0;
  for (unsigned int i=0;i<m_maxn+1;++i)
    m_histos.push_back(new Histogram(type,m_xmin,m_xmax,m_nbins));
}

Jet_Observable_Base::~Jet_Observable_Base()
{
  for (size_t i=0;i<m_histos.size();++i) delete m_histos[i];
}

void Jet_Observable_Base::Evaluate(const Particle_List & pl,double weight, double ncount)
{
  if ((m_mode==1 && pl.size()>=m_minn) ||
      (m_mode==2 && pl.size()==m_minn)) {
    // fill
    size_t i=1;
    m_histos[0]->Insert(0.,0.,ncount);
    for (Particle_List::const_iterator it=pl.begin();it!=pl.end();++it,++i) {
      double value=Calc(*it);
      m_histos[0]->Insert(value,weight,0);
      if (i<=m_maxn) m_histos[i]->Insert(value,weight,ncount);
    }
    for (; i<m_histos.size();++i) { 
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
  else {
    // fill with 0
    m_histos[0]->Insert(0.,0.,ncount);
    for (size_t i=1; i<m_histos.size();++i) {
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
}

void Jet_Observable_Base::Evaluate(const Blob_List & blobs,double value, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}


void Jet_Observable_Base::EvaluateNLOcontrib(double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  if ((m_mode==1 && pl->size()>=m_minn) ||
      (m_mode==2 && pl->size()==m_minn)) {
    // fill
    size_t i=1;
    for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it,++i) {
      double value=Calc(*it);
      m_histos[0]->InsertMCB(value,weight,ncount);
      if (i<=m_maxn) m_histos[i]->InsertMCB(value,weight,ncount);
    }
    for (; i<m_histos.size();++i) { 
      m_histos[i]->InsertMCB(0.,0.,ncount);
    }
  }
  else {
    // fill with 0
    m_histos[0]->InsertMCB(0.,0.,ncount);
    for (size_t i=1; i<m_histos.size();++i) {
      m_histos[i]->InsertMCB(0.,0.,ncount);
    }
  }
}

void Jet_Observable_Base::EvaluateNLOevt()
{
  for (size_t i=0; i<m_histos.size();++i) { 
    m_histos[i]->FinishMCB();
  }
}

void Jet_Observable_Base::EndEvaluation(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Observable_Base::Restore(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
  }
}

void Jet_Observable_Base::Output(const std::string & pname) {
  ATOOLS::MakeDir(pname); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i;
    s<<".dat"; 
    s>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+fname).c_str());
  }
}

Primitive_Observable_Base & Jet_Observable_Base::operator+=(const Primitive_Observable_Base & ob)
{
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    msg_Error()<<"ERROR: in Jet_Observable_Base::operator+=  in"<<m_name<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    return *this;
  }

  Jet_Observable_Base * job = ((Jet_Observable_Base*)(&ob));

  if (m_histos.size()==job->m_histos.size()) {
    for (size_t i=0; i<m_histos.size();++i) {
      (*m_histos[i])+=(*job->m_histos[i]);
    }
  }
  return *this;
}

void Jet_Observable_Base::Reset()
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}



Two_Jet_Observable_Base::Two_Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
						 unsigned int mode,unsigned int minn,unsigned int maxn, 
						 const std::string & lname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), m_mode(mode), m_minn(minn), m_maxn(maxn)
{
  m_listname = lname;
  m_name     = std::string("jet_");
  if (lname!=finalstate_list) m_name=lname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }
  if (p_histo) delete p_histo; 
  p_histo =  0;
  unsigned int num = (m_maxn*m_maxn-m_maxn)/2;
  for (unsigned int i=0;i<num+1;++i)
    m_histos.push_back(new Histogram(type,m_xmin,m_xmax,m_nbins));

  p_minpts = new double[maxn]; p_maxpts = new double[maxn];
  for (unsigned int i=0;i<maxn;i++) { p_minpts[i]=0.; p_maxpts[i]=1.e12; }
}

Two_Jet_Observable_Base::~Two_Jet_Observable_Base()
{
  for (size_t i=0;i<m_histos.size();++i) delete m_histos[i];
  delete [] p_minpts;
  delete [] p_maxpts;
}

void Two_Jet_Observable_Base::Evaluate(const Particle_List & pl,
				       double weight, double ncount)
{
  if((m_mode==1 && pl.size()>=m_minn) || (m_mode==2 && pl.size()==m_minn)) {
    //fill
    size_t i=1, jet1=1, jet2=1;
    size_t j1max=Min(pl.size()-1,size_t(m_maxn)-1), j2max=j1max+1;
    for(Particle_List::const_iterator it1=pl.begin();
	jet1<=j1max; ++it1, ++jet1) {
      jet2=jet1+1;
      for(Particle_List::const_iterator it2=it1+1;
	  jet2<=j2max; ++it2, ++i, ++jet2) {
	double value=Calc(*it1,*it2,jet1-1,jet2-1);
	m_histos[0]->Insert(value,weight,ncount);
	msg_Debugging()<<"2-jet obs '"<<m_name<<"': j#1 = "<<jet1<<", j#2 = "<<jet2
		       <<", idx = "<<i<<"\n";
	m_histos[i]->Insert(value,weight,ncount);
      }
      for(; jet2<=size_t(m_maxn); ++i, ++jet2) {
	m_histos[0]->Insert(0.,0.,ncount);
	m_histos[i]->Insert(0.,0.,ncount);
      }
    }
    for (; i<m_histos.size(); ++i) {
      m_histos[0]->Insert(0.,0.,ncount);
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
  else {
    //fill with 0
    for(size_t i=1; i<m_histos.size(); ++i) {
      m_histos[0]->Insert(0.,0.,ncount);
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
}


void Two_Jet_Observable_Base::Evaluate(const Blob_List & blobs,double value, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}

void Two_Jet_Observable_Base::EvaluateNLOcontrib(double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  if ((m_mode==1 && pl->size()>=m_minn) ||
      (m_mode==2 && pl->size()==m_minn)) {
    //fill
    size_t i=1, jet1=1, jet2=1;
    size_t j1max=Min(pl->size()-1,size_t(m_maxn)-1), j2max=j1max+1;
    for(Particle_List::const_iterator it1=pl->begin();
	jet1<=j1max; ++it1, ++jet1) {
      jet2=jet1+1;
      for(Particle_List::const_iterator it2=it1+1;
	  jet2<=j2max; ++it2, ++i, ++jet2) {
	double value=Calc(*it1,*it2,jet1-1,jet2-1);
	m_histos[0]->InsertMCB(value,weight,ncount);
	msg_Debugging()<<"2-jet obs '"<<m_name<<"': j#1 = "<<jet1<<", j#2 = "<<jet2
		       <<", idx = "<<i<<"\n";
	m_histos[i]->InsertMCB(value,weight,ncount);
      }
      for(; jet2<=size_t(m_maxn); ++i, ++jet2) {
	m_histos[0]->InsertMCB(0.,0.,ncount);
	m_histos[i]->InsertMCB(0.,0.,ncount);
      }
    }
    for (; i<m_histos.size(); ++i) {
      m_histos[0]->InsertMCB(0.,0.,ncount);
      m_histos[i]->InsertMCB(0.,0.,ncount);
    }
  }
  else {
    // fill with 0
    for (size_t i=0; i<m_histos.size();++i) {
      m_histos[0]->InsertMCB(0.,0.,ncount);
      m_histos[i]->InsertMCB(0.,0.,ncount);
    }
  }
}

void Two_Jet_Observable_Base::EvaluateNLOevt()
{
  if (m_histos.size()<2) return;
  for (size_t i=0; i<m_histos.size();++i) { 
    m_histos[i]->FinishMCB();
  }
}

void Two_Jet_Observable_Base::EndEvaluation(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Two_Jet_Observable_Base::Restore(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
  }
}

void Two_Jet_Observable_Base::Output(const std::string & pname) {
  ATOOLS::MakeDir(pname); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i;
    s<<".dat"; 
    s>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+fname).c_str());
  }
}

Primitive_Observable_Base & Two_Jet_Observable_Base::operator+=(const Primitive_Observable_Base & ob)
{
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    std::cout<<" ERROR: in Two_Jet_Observable_Base::operator+=  in"<<m_name<<std::endl;
    return *this;
  }

  Two_Jet_Observable_Base * jdrd = ((Two_Jet_Observable_Base*)(&ob));

  if (m_histos.size()==jdrd->m_histos.size()) {
    for (size_t i=0; i<m_histos.size();++i) {
      (*m_histos[i])+=(*jdrd->m_histos[i]);
    }
  }
  return *this;
}

void Two_Jet_Observable_Base::Reset()
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}

void Two_Jet_Observable_Base::SetPTRange(const unsigned int jetno,const double minpt,const double maxpt)
{
  if (!(/*jetno>=m_minn && */ jetno<=m_maxn)) {
    msg_Error()<<"Potential Error in Two_Jet_Observable_Base::SetMinPT("<<jetno<<")"<<std::endl
	       <<"   Out of bounds : "<<m_minn<<" ... "<<m_maxn<<", will continue."<<std::endl;
    return;
  }
  p_minpts[jetno-1] = minpt; 
  p_maxpts[jetno-1] = maxpt; 
}




//########################################################################################
//########################################################################################
//########################################################################################

DEFINE_OBSERVABLE_GETTER(Jet_Rapidity_Distribution,Jet_Rapidity_Distribution_Getter,"JetRap")

Jet_Rapidity_Distribution::Jet_Rapidity_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					   unsigned int mode,unsigned int minn,unsigned int maxn, 
					   const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="y_";
}


double Jet_Rapidity_Distribution::Calc(const Particle * p)
{
  Vec4D mom=p->Momentum();
  return mom.Y();
}

Primitive_Observable_Base * Jet_Rapidity_Distribution::Copy() const 
{
  return new Jet_Rapidity_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

//########################################################################################
//########################################################################################
//########################################################################################


DEFINE_OBSERVABLE_GETTER(Jet_Eta_Distribution,Jet_Eta_Distribution_Getter,"JetEta")

Jet_Eta_Distribution::Jet_Eta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					   unsigned int mode,unsigned int minn,unsigned int maxn, 
					   const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="eta_";
}


double Jet_Eta_Distribution::Calc(const Particle * p)
{
  Vec4D mom=p->Momentum();

  double pt2=sqr(mom[1])+sqr(mom[2]);
  double pp =sqrt(pt2+sqr(mom[3]));
  double pz =dabs(mom[3]);
  double sn =mom[3]/pz;
  if (pt2<1.e-10*pp*pp) {
    return sn*20.;
  }
  return sn*0.5*log(sqr(pp+pz)/pt2);
}

Primitive_Observable_Base * Jet_Eta_Distribution::Copy() const 
{
  return new Jet_Eta_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Jet_Phi_Distribution,Jet_Phi_Distribution_Getter,"JetPhi")

Jet_Phi_Distribution::Jet_Phi_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					   unsigned int mode,unsigned int minn,unsigned int maxn, 
					   const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="phi_";
}

Primitive_Observable_Base * Jet_Phi_Distribution::Copy() const 
{
  return new Jet_Phi_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

double Jet_Phi_Distribution::Calc(const Particle * p)
{
  Vec4D mom=p->Momentum();
  double phi = mom.Phi();
  return phi;
}

DEFINE_OBSERVABLE_GETTER(Jet_PT_Distribution,Jet_PT_Distribution_Getter,"JetPT")

Jet_PT_Distribution::Jet_PT_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="pt_";
}


double Jet_PT_Distribution::Calc(const Particle * p)
{
  Vec4D mom=p->Momentum();
//   PRINT_INFO(sqrt(sqr(mom[1])+sqr(mom[2])));
  return sqrt(sqr(mom[1])+sqr(mom[2]));
}

Primitive_Observable_Base * Jet_PT_Distribution::Copy() const 
{
  return new Jet_PT_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Jet_IPT2_Distribution,Jet_IPT2_Distribution_Getter,"JetIPT2")

Jet_IPT2_Distribution::Jet_IPT2_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="ipt2_";
}


double Jet_IPT2_Distribution::Calc(const Particle * p)
{
  return 1.0/Max(1.0e-12,p->Momentum().PPerp2());
}

Primitive_Observable_Base * Jet_IPT2_Distribution::Copy() const 
{
  return new Jet_IPT2_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Jet_ET_Distribution,Jet_ET_Distribution_Getter,"JetET")

Jet_ET_Distribution::Jet_ET_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="Et_";
}


double Jet_ET_Distribution::Calc(const Particle * p)
{
  Vec4D mom=p->Momentum();
  double pt2 = sqr(mom[1])+sqr(mom[2]);
  double p2  = sqr(mom[3])+pt2;
  return mom[0]*sqrt(pt2/p2);
}

Primitive_Observable_Base * Jet_ET_Distribution::Copy() const 
{
  return new Jet_ET_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Jet_E_Distribution,Jet_E_Distribution_Getter,"JetE")

Jet_E_Distribution::Jet_E_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name+="E_";
}


double Jet_E_Distribution::Calc(const Particle * p)
{
  Vec4D mom=p->Momentum();
  return mom[0];
}

Primitive_Observable_Base * Jet_E_Distribution::Copy() const 
{
  return new Jet_E_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

DECLARE_GETTER(Jet_Differential_Rates_Getter,"JetDRate",
	       Primitive_Observable_Base,Argument_Matrix);	

DEFINE_GETTER_METHOD(Jet_Differential_Rates,Jet_Differential_Rates_Getter)

void Jet_Differential_Rates_Getter::PrintInfo(std::ostream &str,const size_t width) const	
{ 
  str<<"min max bins mode nmin nmax Lin|LinErr|Log|LogErr [list] -> Finder 93 .."; 
}


Jet_Differential_Rates::Jet_Differential_Rates(unsigned int type,double xmin,double xmax,int nbins,
					       unsigned int mode,unsigned int minn,unsigned int maxn, 
					       const std::string & listname,
					       const std::string & reflistname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
//   std::cout<<"Jet_Differential_Rates("<<type<<", "<<xmin<<", "<<xmax<<", "<<nbins<<",\n"
// 	   <<mode<<", "<<minn<<", "<<maxn<<",\n"
// 	   <<listname<<", "<<reflistname<<")\n";
  if (reflistname=="") {
    m_reflistname = listname;
    m_name=listname+"_KtJetrates(1)jet_";
  }
  else {
    m_reflistname = reflistname;
    m_name=listname+"_"+reflistname+"_KtJetrates(1)jet_";
  } 
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }
}


void Jet_Differential_Rates::Evaluate(const Blob_List & blobs,double weight, double ncount)
{
  Blob_Data_Base * ktdrs=(*p_ana)["KtDeltaRs"];
  std::string key="KtJetrates(1)"+m_listname;
  if (ktdrs) {
//     std::vector<double> * drs=ktdrs->Get<std::vector<double> *>();
    /*
    MyStrStream str;
    str<<"KtJetrates("<<(*drs)[0]<<")"<<m_listname;
    str>>key;
    */
    key="KtJetrates(1)"+m_listname;
  }

  Blob_Data_Base * rates=(*p_ana)[key];
  if (!rates) {
    msg_Out()<<"WARNING in Jet_Differential_Rates::Evaluate : "<<key<<" not found "<<std::endl;
    return;
  }
  Particle_List * pl=p_ana->GetParticleList(m_reflistname);
  if (!pl) {
    msg_Out()<<"WARNING in Jet_Differential_Rates::Evaluate : "<<m_reflistname<<" not found "<<std::endl;
    return;
  }

  std::vector<double> * jd=rates->Get<std::vector<double> *>();

  size_t j=jd->size();
  // plot only selected events 
  if (pl->size()==0) j=0;

  for (size_t i=0; i<m_histos.size();++i) {
    if (j>0) {
      --j;
      m_histos[i]->Insert(sqrt((*jd)[j]),weight,ncount);
    }
    else {
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
}

void Jet_Differential_Rates::EvaluateNLOcontrib(double weight, double ncount)
{
  Blob_Data_Base * ktdrs=(*p_ana)["KtDeltaRs"];
  std::string key="KtJetrates(1)"+m_listname;
  if (ktdrs) {
//     std::vector<double> * drs=ktdrs->Get<std::vector<double> *>();
    /*
    MyStrStream str;
    str<<"KtJetrates("<<(*drs)[0]<<")"<<m_listname;
    str>>key;
    */
    key="KtJetrates(1)"+m_listname;
  }

  Blob_Data_Base * rates=(*p_ana)[key];
  if (!rates) {
    msg_Out()<<"WARNING in Jet_Differential_Rates::Evaluate : "<<key<<" not found "<<std::endl;
    return;
  }
  Particle_List * pl=p_ana->GetParticleList(m_reflistname);
  if (!pl) {
    msg_Out()<<"WARNING in Jet_Differential_Rates::Evaluate : "<<m_reflistname<<" not found "<<std::endl;
    return;
  }

  std::vector<double> * jd=rates->Get<std::vector<double> *>();

  size_t j=jd->size();
  // plot only selected events 
  if (pl->size()==0) j=0;

  for (size_t i=0; i<m_histos.size();++i) {
    if (j>0) {
      --j;
      m_histos[i]->InsertMCB(sqrt((*jd)[j]),weight,ncount);
    }
    else {
      m_histos[i]->InsertMCB(0.,weight,ncount);
    }
  }
}

double Jet_Differential_Rates::Calc(const Particle *) 
{
  return 0.;
}


Primitive_Observable_Base * Jet_Differential_Rates::Copy() const 
{
  if (m_listname==m_reflistname)
    return new Jet_Differential_Rates(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,
				      m_listname);
  else
    return new Jet_Differential_Rates(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,
				      m_listname,m_reflistname);
}




//##############################################################################
//##############################################################################
//##############################################################################


DEFINE_OBSERVABLE_GETTER(Jet_DeltaR_Distribution,
			 Jet_DeltaR_Distribution_Getter,"JetDR")

////////////////////////////////////////////////////////////////////////////////

Jet_DeltaR_Distribution::Jet_DeltaR_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						 unsigned int mode,unsigned int minn,unsigned int maxn, 
						 const std::string & lname) :
  Two_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="dR2_";
}

Primitive_Observable_Base * Jet_DeltaR_Distribution::Copy() const 
{
  Jet_DeltaR_Distribution * jdr =
    new Jet_DeltaR_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  for (unsigned int i=0;i<m_maxn;i++) jdr->SetPTRange(i+1,p_minpts[i],p_maxpts[i]);
  return jdr;
}

double Jet_DeltaR_Distribution::Calc(const Particle * p1,const Particle * p2,
				     const int jet1,const int jet2)
{
  Vec4D mom1=p1->Momentum();
  Vec4D mom2=p2->Momentum();
  
  double pt1  = mom1.PPerp();
  double pt2  = mom2.PPerp();
  if (pt1<p_minpts[jet1] || pt2<p_minpts[jet2] ||
      pt1>p_maxpts[jet1] || pt2>p_maxpts[jet2]) return 0.;
  double dphi = acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
  double deta = mom1.Eta()-mom2.Eta();
  return sqrt(sqr(deta) + sqr(dphi)); 
}

//----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(Jet_DeltaEta_Distribution,
			 Jet_DeltaEta_Distribution_Getter,"JetDEta")

Jet_DeltaEta_Distribution::Jet_DeltaEta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						     unsigned int mode,unsigned int minn,unsigned int maxn, 
						     const std::string & lname) :
  Two_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="deta2_";
}

Primitive_Observable_Base * Jet_DeltaEta_Distribution::Copy() const 
{
  Jet_DeltaEta_Distribution * jde =
    new Jet_DeltaEta_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  for (unsigned int i=0;i<m_maxn;i++) jde->SetPTRange(i+1,p_minpts[i],p_maxpts[i]);
  return jde;
}

double Jet_DeltaEta_Distribution::Calc(const Particle * p1,const Particle * p2,
				       const int jet1,const int jet2)
{
  Vec4D mom1 = p1->Momentum();
  Vec4D mom2 = p2->Momentum();
  double pt1 = mom1.PPerp(), pt2 = mom2.PPerp();
  if (pt1<p_minpts[jet1] || pt2<p_minpts[jet2] ||
      pt1>p_maxpts[jet1] || pt2>p_maxpts[jet2]) return 0.;
  
  return dabs((mom1.Eta()-mom2.Eta()));
}
//----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(Jet_DeltaY_Distribution,
			 Jet_DeltaY_Distribution_Getter,"JetDY")

Jet_DeltaY_Distribution::Jet_DeltaY_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						     unsigned int mode,unsigned int minn,unsigned int maxn, 
						     const std::string & lname) :
  Two_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="dy2_";
}

Primitive_Observable_Base * Jet_DeltaY_Distribution::Copy() const 
{
  Jet_DeltaY_Distribution * jde =
    new Jet_DeltaY_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  return jde;
}

double Jet_DeltaY_Distribution::Calc(const Particle * p1,const Particle * p2,
				       const int jet1,const int jet2)
{
  Vec4D mom1 = p1->Momentum();
  Vec4D mom2 = p2->Momentum();
  
  double y1(mom1.Y()), y2(mom2.Y());
  if (!(y1>=0.0) && !(y1<0.0)) return 0.0;
  if (!(y2>=0.0) && !(y2<0.0)) return 0.0;
  return dabs(y1-y2);
}
//----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(Jet_DeltaPhi_Distribution,
			 Jet_DeltaPhi_Distribution_Getter,"JetDPhi")

Jet_DeltaPhi_Distribution::Jet_DeltaPhi_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						     unsigned int mode,unsigned int minn,unsigned int maxn, 
						     const std::string & lname) :
  Two_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="dphi2_";
}

Primitive_Observable_Base * Jet_DeltaPhi_Distribution::Copy() const 
{
  Jet_DeltaPhi_Distribution * jdp =
    new Jet_DeltaPhi_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  for (unsigned int i=0;i<m_maxn;i++) jdp->SetPTRange(i+1,p_minpts[i],p_maxpts[i]);
  return jdp;
}

double Jet_DeltaPhi_Distribution::Calc(const Particle * p1,const Particle * p2,
				       const int jet1,const int jet2)
{
  Vec4D mom1=p1->Momentum();
  Vec4D mom2=p2->Momentum();
  double pt1 = mom1.PPerp(), pt2 = mom2.PPerp();
  if (pt1<p_minpts[jet1] || pt2<p_minpts[jet2] ||
      pt1>p_maxpts[jet1] || pt2>p_maxpts[jet2]) return 0.;
  
  return acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
}
//----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(Jet_DiMass_Distribution,
			 Jet_DiMass_Distribution_Getter,"JetDiMass")

Jet_DiMass_Distribution::Jet_DiMass_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						     unsigned int mode,unsigned int minn,unsigned int maxn, 
						     const std::string & lname) :
  Two_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="DJM";
}

Primitive_Observable_Base * Jet_DiMass_Distribution::Copy() const 
{
  Jet_DiMass_Distribution * jdp =
    new Jet_DiMass_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  for (unsigned int i=0;i<m_maxn;i++) jdp->SetPTRange(i+1,p_minpts[i],p_maxpts[i]);
  return jdp;
}

double Jet_DiMass_Distribution::Calc(const Particle * p1,const Particle * p2,
				       const int jet1,const int jet2)
{
  Vec4D mom1=p1->Momentum();
  Vec4D mom2=p2->Momentum();
  double pt1 = mom1.PPerp(), pt2 = mom2.PPerp();
  if (pt1<p_minpts[jet1] || pt2<p_minpts[jet2] ||
      pt1>p_maxpts[jet1] || pt2>p_maxpts[jet2]) return 0.;
  
  return (mom1+mom2).Abs();
}

Three_Jet_Observable_Base::Three_Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
                                                 unsigned int mode,unsigned int minn,unsigned int maxn, 
                                                 const std::string & lname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), m_mode(mode), m_minn(minn), m_maxn(maxn)
{
  m_listname = lname;
  m_name     = std::string("jet_");
  if (lname!=finalstate_list) m_name=lname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }
  if (p_histo) delete p_histo;
  p_histo =  0;
  unsigned int num = m_maxn*(m_maxn-1)*(m_maxn-2)/6;
  for (unsigned int i=0;i<num+1;++i)
    m_histos.push_back(new Histogram(type,m_xmin,m_xmax,m_nbins));

}

Three_Jet_Observable_Base::~Three_Jet_Observable_Base()
{
  for (size_t i=0;i<m_histos.size();++i) delete m_histos[i];
}

void Three_Jet_Observable_Base::Evaluate(const Particle_List & pl,double weight, double ncount)
{
  if ((m_mode==1 && pl.size()>=m_minn) ||
      (m_mode==2 && pl.size()==m_minn)) {
    // fill
    size_t i=1;
    int jet1=0;
    for (Particle_List::const_iterator it1=pl.begin();it1!=pl.end();++it1,++jet1) {
      int jet2=jet1+1;
      for (Particle_List::const_iterator it2=it1+1;it2!=pl.end() && i<=m_maxn*(m_maxn-1)*(m_maxn-2)/6;++it2,++i,++jet2) {
        int jet3=jet2+1;
        for (Particle_List::const_iterator it3=it2+1;it3!=pl.end() && i<=m_maxn*(m_maxn-1)*(m_maxn-2)/6;++it3,++i,++jet3) {
          double value=Calc(*it1,*it2,*it3,jet1,jet2,jet3);
          m_histos[0]->Insert(value,weight,ncount);
          m_histos[i]->Insert(value,weight,ncount);
        }
      }
    }
    for (; i<m_histos.size();++i) { 
      m_histos[0]->Insert(0.,0.,ncount);
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
  else {
    // fill with 0
    for (size_t i=0; i<m_histos.size();++i) {
      m_histos[0]->Insert(0.,0.,ncount);
      m_histos[i]->Insert(0.,0.,ncount);
    }
  }
}

void Three_Jet_Observable_Base::Evaluate(const Blob_List & blobs,double value, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}

void Three_Jet_Observable_Base::EndEvaluation(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Three_Jet_Observable_Base::Restore(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
  }
}

void Three_Jet_Observable_Base::Output(const std::string & pname) {
  ATOOLS::MakeDir(pname); 
  for (size_t i=0; i<m_histos.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i;
    s<<".dat"; 
    s>>fname;
    m_histos[i]->Output((pname+std::string("/")+m_name+fname).c_str());
  }
}

Primitive_Observable_Base & Three_Jet_Observable_Base::operator+=(const Primitive_Observable_Base & ob)
{
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    std::cout<<" ERROR: in Three_Jet_Observable_Base::operator+=  in"<<m_name<<std::endl;
    return *this;
  }

  Three_Jet_Observable_Base * jdrd = ((Three_Jet_Observable_Base*)(&ob));

  if (m_histos.size()==jdrd->m_histos.size()) {
    for (size_t i=0; i<m_histos.size();++i) {
      (*m_histos[i])+=(*jdrd->m_histos[i]);
    }
  }
  return *this;
}

void Three_Jet_Observable_Base::Reset()
{
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Reset();
  }  
}

DEFINE_OBSERVABLE_GETTER(Eta_3_Prime,
			 Eta_3_Prime_Getter,"Eta3Prime")

Eta_3_Prime::Eta_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & lname) :
  Three_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="eta3p_";
}

Primitive_Observable_Base * Eta_3_Prime::Copy() const 
{
  Eta_3_Prime * jde =
    new Eta_3_Prime(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  return jde;
}

double Eta_3_Prime::Calc(const Particle * p1,const Particle * p2,
			 const Particle * p3,const int jet1,const int jet2,
			 const int jet3)
{
  Vec4D mom1 = p1->Momentum();
  Vec4D mom2 = p2->Momentum();
  Vec4D mom3 = p3->Momentum();
//   PRINT_INFO(mom1<<mom2<<mom3);
  
  double y1(mom1.Eta()), y2(mom2.Eta()), y3(mom3.Eta());
  if (!(y1>=0.0) && !(y1<0.0)) return 0.0;
  if (!(y2>=0.0) && !(y2<0.0)) return 0.0;
  if (!(y3>=0.0) && !(y3<0.0)) return 0.0;
  return y3-(y1+y2)/2.0;
}

DEFINE_OBSERVABLE_GETTER(Y_3_Prime,
			 Y_3_Prime_Getter,"Y3Prime")

Y_3_Prime::Y_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & lname) :
  Three_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="y3p_";
}

Primitive_Observable_Base * Y_3_Prime::Copy() const 
{
  Y_3_Prime * jde =
    new Y_3_Prime(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  return jde;
}

double Y_3_Prime::Calc(const Particle * p1,const Particle * p2,
			 const Particle * p3,const int jet1,const int jet2,
			 const int jet3)
{
  Vec4D mom1 = p1->Momentum();
  Vec4D mom2 = p2->Momentum();
  Vec4D mom3 = p3->Momentum();
  
  double y1(mom1.Y()), y2(mom2.Y()), y3(mom3.Y());
  if (!(y1>=0.0) && !(y1<0.0)) return 0.0;
  if (!(y2>=0.0) && !(y2<0.0)) return 0.0;
  if (!(y3>=0.0) && !(y3<0.0)) return 0.0;
  return y3-(y1+y2)/2.0;
}

DEFINE_OBSERVABLE_GETTER(Phi_3_Prime,
			 Phi_3_Prime_Getter,"Phi3Prime")

Phi_3_Prime::Phi_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & lname) :
  Three_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="phi3p_";
}

Primitive_Observable_Base * Phi_3_Prime::Copy() const 
{
  Phi_3_Prime * jde =
    new Phi_3_Prime(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  return jde;
}

double Phi_3_Prime::Calc(const Particle * p1,const Particle * p2,
			 const Particle * p3,const int jet1,const int jet2,
			 const int jet3)
{
  Vec4D mom1 = p1->Momentum();
  Vec4D mom2 = p2->Momentum();
  Vec4D mom3 = p3->Momentum();
 
  Vec3D n1(mom1), n2(mom2), n3(mom3);
  n1=1.0/n1.Abs()*n1;
  n2=1.0/n2.Abs()*n2;
  n3=1.0/n3.Abs()*n3;

  double phi(acos(n3*cross(n1,n2)));
 
  if (!(phi>=0.0) && !(phi<0.0)) return 0.0;
  return phi;
}

DEFINE_OBSERVABLE_GETTER(CosPhi_3_Prime,
			 CosPhi_3_Prime_Getter,"CosPhi3Prime")

CosPhi_3_Prime::CosPhi_3_Prime(unsigned int type,double xmin,double xmax,int nbins,
			 unsigned int mode,unsigned int minn,unsigned int maxn, 
			 const std::string & lname) :
  Three_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname)
{
  m_name+="cosphi3p_";
}

Primitive_Observable_Base * CosPhi_3_Prime::Copy() const 
{
  CosPhi_3_Prime * jde =
    new CosPhi_3_Prime(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  return jde;
}

double CosPhi_3_Prime::Calc(const Particle * p1,const Particle * p2,
			 const Particle * p3,const int jet1,const int jet2,
			 const int jet3)
{
  Vec4D mom1 = p1->Momentum();
  Vec4D mom2 = p2->Momentum();
  Vec4D mom3 = p3->Momentum();
 
  Vec3D n1(mom1), n2(mom2), n3(mom3);
  n1=1.0/n1.Abs()*n1;
  n2=1.0/n2.Abs()*n2;
  n3=1.0/n3.Abs()*n3;

  double cosphi(n3*cross(n1,n2));
 
  if (!(cosphi>=0.0) && !(cosphi<0.0)) return 0.0;
  return cosphi;
}

DEFINE_OBSERVABLE_GETTER(Jet_Alpha_Distribution,
                         Jet_Alpha_Distribution_Getter,"JetAlpha")

Jet_Alpha_Distribution::Jet_Alpha_Distribution(unsigned int type,double
xmin,
                                               double xmax,int nbins,
                                               unsigned int mode,
                                               unsigned int minn,
                                               unsigned int maxn,
                                               const std::string & lname)
  : Two_Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname) {
  m_name+="Alpha_";
}

Primitive_Observable_Base * Jet_Alpha_Distribution::Copy() const
{
  Jet_Alpha_Distribution * ja =
    new 
    Jet_Alpha_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
  for (unsigned int i=0;i<m_maxn;i++) 
    ja->SetPTRange(i+1,p_minpts[i],p_maxpts[i]);
  return ja;
}

double Jet_Alpha_Distribution::Calc(const Particle* p1, const Particle* p2,
                                    const int jet1, const int jet2) {
  Vec4D mom1=p1->Momentum();
  Vec4D mom2=p2->Momentum(); 

  double dH = mom1.Eta()/fabs(mom1.Eta())*(mom2.Eta()-mom1.Eta());

  double pt1  = mom1.PPerp();
  double pt2  = mom2.PPerp();
  if (pt1<p_minpts[jet1] || pt2<p_minpts[jet2] ||
      pt1>p_maxpts[jet1] || pt2>p_maxpts[jet2]) return 0.;
  double dphi = acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));

  return atan(dH/dphi)/M_PI*180.;
}
