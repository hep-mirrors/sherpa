#include "Jet_Observables.H"
#include "MyStrStream.H"
#include "Primitive_Analysis.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Jet_Observable_Base::Jet_Observable_Base(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & lname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL), m_mode(mode), m_minn(minn), m_maxn(maxn)
{
  m_listname=lname;
  m_name  = std::string("jet_");
  if (lname!="Analysed") m_name=lname+std::string("_")+m_name;
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_minn<<"_"<<m_mode<<"_";
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
    // fuelle
    size_t i=1;
    for (Particle_List::const_iterator it=pl.begin();it!=pl.end() && i<m_maxn;++it,++i) {
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
    // fuelle mit 0
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
  mkdir((pname).c_str(),mode_dir); 
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

Jet_Eta_Distribution::Jet_Eta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					   unsigned int mode,unsigned int minn,unsigned int maxn, 
					   const std::string & lname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname) 
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

Jet_PT_Distribution::Jet_PT_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & lname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname) 
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


Jet_E_Distribution::Jet_E_Distribution(unsigned int type,double xmin,double xmax,int nbins,
					 unsigned int mode,unsigned int minn,unsigned int maxn, 
					 const std::string & lname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname) 
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


Jet_Differential_Rates::Jet_Differential_Rates(unsigned int type,double xmin,double xmax,int nbins,
					       unsigned int mode,unsigned int minn,unsigned int maxn, 
					       const std::string & lname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,lname) 
{
  m_name="KtJetrates(1)"+m_name;
}


void Jet_Differential_Rates::Evaluate(const Blob_List & blobs,double weight, int ncount)
{
  Blob_Data_Base * ktdrs=(*p_ana)["KtDeltaRs"];
  std::string key="KtJetrates(1)"+m_listname;
  if (ktdrs) {
    std::vector<double> * drs=ktdrs->Get<std::vector<double> *>();
    MyStrStream str;
    str<<"KtJetrates("<<(*drs)[0]<<")"<<m_listname;
    str>>key;
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
