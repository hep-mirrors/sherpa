#include "Jet_Observables.H"
#include "MyStrStream.H"
#include "Primitive_Analysis.H"
#include "Shell_Tools.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class,class Getter>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    std::string list=parameters[0].size()>7?parameters[0][7]:"Analysed";
    return new Class(10*(int)(parameters[0][6]=="Log"),
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
  std::string list="Analysed", scale="Lin";
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
  return new Class((scale=="Log")*10,min,max,bins,mode,nmin,nmax,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME,TAG)				\
  Primitive_Observable_Base *const					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS,NAME>(parameters); }

#define DEFINE_PRINT_METHOD(CLASS,NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins mode nmin nmax Lin|Log [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME,TAG);					\
  DEFINE_PRINT_METHOD(CLASS,NAME)

using namespace ATOOLS;

Jet_Observable_Base::Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL), m_mode(mode), m_minn(minn), m_maxn(maxn)
{
  m_listname=listname;
  m_name  = std::string("jet_");
  if (listname!="Analysed") m_name=listname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }

  p_histo =  0;
  for (unsigned int i=0;i<m_maxn+1;++i)
    m_histos.push_back(new Histogram(type,m_xmin,m_xmax,m_nbins));
}

void Jet_Observable_Base::Evaluate(const Particle_List & pl,double weight, int ncount)
{
  if ((m_mode==1 && pl.size()>=m_minn) ||
      (m_mode==2 && pl.size()==m_minn)) {
    // fill
    size_t i=1;
    for (Particle_List::const_iterator it=pl.begin();it!=pl.end() && i<=m_maxn;++it,++i) {
      double value=Calc(*it);
      m_histos[0]->Insert(value,weight,ncount);
      m_histos[i]->Insert(value,weight,ncount);
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

void Jet_Observable_Base::Evaluate(const Blob_List & blobs,double value, int ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}

void Jet_Observable_Base::EndEvaluation(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Jet_Observable_Base::Output(const std::string & pname) {
  int  mode_dir = 448;
  ATOOLS::MakeDir((pname).c_str(),mode_dir); 
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
    msg.Error()<<"ERROR: in Jet_Observable_Base::operator+=  in"<<m_name<<std::endl
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
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL), m_mode(mode), m_minn(minn), m_maxn(maxn)
{
  m_listname = lname;
  m_name     = std::string("jet_");
  if (lname!="Analysed") m_name=lname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }

  p_histo =  0;
  unsigned int num = (m_maxn*m_maxn-m_maxn)/2;
  for (unsigned int i=0;i<num+1;++i)
    m_histos.push_back(new Histogram(type,m_xmin,m_xmax,m_nbins));

  p_minpts = new double[maxn]; p_maxpts = new double[maxn];
  for (unsigned int i=0;i<maxn;i++) { p_minpts[i]=0.; p_maxpts[i]=1.e12; }
}


void Two_Jet_Observable_Base::Evaluate(const Particle_List & pl,double weight, int ncount)
{
  if ((m_mode==1 && pl.size()>=m_minn) ||
      (m_mode==2 && pl.size()==m_minn)) {
    // fill
    size_t i=1;
    int jet1=0,jet2=0;
    for (Particle_List::const_iterator it1=pl.begin();it1!=pl.end();++it1,++jet1) {
      for (Particle_List::const_iterator it2=it1+1;it2!=pl.end() && i<=(sqr(m_maxn)-m_maxn)/2;++it2,++i,++jet2) {
	double value=Calc(*it1,*it2,jet1,jet2);
	m_histos[0]->Insert(value,weight,ncount);
	m_histos[i]->Insert(value,weight,ncount);
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

void Two_Jet_Observable_Base::Evaluate(const Blob_List & blobs,double value, int ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}

void Two_Jet_Observable_Base::EndEvaluation(double scale) {
  for (size_t i=0; i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    if (scale!=1.) m_histos[i]->Scale(scale);
    m_histos[i]->Output();
  }
}

void Two_Jet_Observable_Base::Output(const std::string & pname) {
  int  mode_dir = 448;
  ATOOLS::MakeDir((pname).c_str(),mode_dir); 
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
  if (!(jetno>=m_minn && jetno<=m_maxn)) {
    msg.Error()<<"Potential Error in Two_Jet_Observable_Base::SetMinPT("<<jetno<<")"<<std::endl
	       <<"   Out of bounds : "<<m_minn<<" ... "<<m_maxn<<", will continue."<<std::endl;
    return;
  }
  p_minpts[jetno-1] = minpt; 
  p_maxpts[jetno-1] = maxpt; 
}




//########################################################################################
//########################################################################################
//########################################################################################

DEFINE_OBSERVABLE_GETTER(Jet_Eta_Distribution,Jet_Eta_Distribution_Getter,"JetEta");

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

DEFINE_OBSERVABLE_GETTER(Jet_PT_Distribution,Jet_PT_Distribution_Getter,"JetPT");

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
  return sqrt(sqr(mom[1])+sqr(mom[2]));
}

Primitive_Observable_Base * Jet_PT_Distribution::Copy() const 
{
  return new Jet_PT_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Jet_E_Distribution,Jet_E_Distribution_Getter,"JetE");

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
	       Primitive_Observable_Base,String_Matrix);	

DEFINE_GETTER_METHOD(Jet_Differential_Rates,Jet_Differential_Rates_Getter,"JetDRate");					

void Jet_Differential_Rates_Getter::PrintInfo(std::ostream &str,const size_t width) const	
{ 
  str<<"min max bins mode nmin nmax Lin|Log [list] -> Finder 93 .."; 
}


Jet_Differential_Rates::Jet_Differential_Rates(unsigned int type,double xmin,double xmax,int nbins,
					       unsigned int mode,unsigned int minn,unsigned int maxn, 
					       const std::string & listname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  m_name="KtJetrates(1)"+m_name;
}


void Jet_Differential_Rates::Evaluate(const Blob_List & blobs,double weight, int ncount)
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
    msg.Out()<<"WARNING in Jet_Differential_Rates::Evaluate : "<<key<<" not found "<<std::endl;
    return;
  }
  std::vector<double> * jd=rates->Get<std::vector<double> *>();

  size_t j=jd->size();
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

double Jet_Differential_Rates::Calc(const Particle *) 
{
  return 0.;
}


Primitive_Observable_Base * Jet_Differential_Rates::Copy() const 
{
  return new Jet_Differential_Rates(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname);
}




//##############################################################################
//##############################################################################
//##############################################################################


DEFINE_OBSERVABLE_GETTER(Jet_DeltaR_Distribution,
			 Jet_DeltaR_Distribution_Getter,"JetDR");

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
			 Jet_DeltaEta_Distribution_Getter,"JetDEta");

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

DEFINE_OBSERVABLE_GETTER(Jet_DeltaPhi_Distribution,
			 Jet_DeltaPhi_Distribution_Getter,"JetDPhi");

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
