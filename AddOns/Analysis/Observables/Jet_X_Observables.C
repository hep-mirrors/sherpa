#include "AddOns/Analysis/Observables/Jet_Observables.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    std::string list=parameters[0].size()>7?parameters[0][7]:finalstate_list;
    std::string reflist=parameters[0].size()>8?parameters[0][8]:"";
    return new Class(HistogramType(parameters[0][6]),
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
  return new Class(HistogramType(scale),min,max,bins,mode,nmin,nmax,list,reflist);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins mode nmin nmax Lin|LinErr|Log|LogErr list reflist"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;

  class Jet_X_Observables : public Jet_Observable_Base {
  protected:
    std::string m_reflistname;
  public:
    Jet_X_Observables(unsigned int type,double xmin,double xmax,int nbins,
			   unsigned int mode, unsigned int minn,unsigned int maxn, 
			   const std::string & =std::string("FinalState"),
			   const std::string & =std::string(""));
    void Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    double Calc(const Particle * p);
    virtual double Calc(const Vec4D &m1,const Vec4D &m2) =0;
  };


Jet_X_Observables::Jet_X_Observables(unsigned int type,double xmin,double xmax,int nbins,
					       unsigned int mode,unsigned int minn,unsigned int maxn, 
					       const std::string & listname,
					       const std::string & reflistname) :
  Jet_Observable_Base(type,xmin,xmax,nbins,mode,minn,maxn,listname) 
{
  if (reflistname=="") {
    m_reflistname = listname;
    m_name=listname+"_";
  }
  else {
    m_reflistname = reflistname;
    m_name=listname+"_"+reflistname+"_";
  } 
  if (m_minn!=0) {
    MyStrStream str;
    str<<m_name<<m_mode<<"_"<<m_minn<<"_";
    str>>m_name;
  }
}

double Jet_X_Observables::Calc(const Particle * p)
{
  return 0.;
}

void Jet_X_Observables::Evaluate(const Blob_List & blobs,double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  if ((m_mode==1 && pl->size()>=m_minn) ||
      (m_mode==2 && pl->size()==m_minn)) {
    Particle_List * pl2=p_ana->GetParticleList(m_reflistname);
    // fill
    size_t i=1;
    m_histos[0]->Insert(0.,0.,ncount);
    Vec4D rmom(0.,0.,0.,0.);
    for (Particle_List::const_iterator it=pl2->begin();it!=pl2->end();++it) {
      rmom+=(*it)->Momentum();
    }

    for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it,++i) {
      double value=Calc((*it)->Momentum(),rmom);
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

void Jet_X_Observables::EvaluateNLOcontrib(double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  if ((m_mode==1 && pl->size()>=m_minn) ||
      (m_mode==2 && pl->size()==m_minn)) {
    Particle_List * pl2=p_ana->GetParticleList(m_reflistname);
    // fill
    size_t i=1;
    m_histos[0]->InsertMCB(0.,0.,ncount);
    Vec4D rmom(0.,0.,0.,0.);
    for (Particle_List::const_iterator it=pl2->begin();it!=pl2->end();++it)
      rmom+=(*it)->Momentum();

    for (Particle_List::const_iterator it=pl->begin();it!=pl->end();++it,++i) {
      double value=Calc((*it)->Momentum(),rmom);
      m_histos[0]->InsertMCB(value,weight,0);
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


////////////////////////////////////////////////////////////////////////////////

  class JetX_DeltaR_Distribution : public Jet_X_Observables {
  protected:
  public:
    JetX_DeltaR_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			      unsigned int mode, unsigned int minn, unsigned int maxn, 
			      const std::string & =std::string("Jets"),
			      const std::string & =std::string("FinalState"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1,const Vec4D &m2);
  };

DEFINE_OBSERVABLE_GETTER(JetX_DeltaR_Distribution,
			 JetX_DeltaR_Distribution_Getter,"JetXDR")

JetX_DeltaR_Distribution::JetX_DeltaR_Distribution(unsigned int type,double xmin,double xmax,int nbins,
						    unsigned int mode,unsigned int minn,unsigned int maxn, 
						    const std::string & lname,const std::string & rname) :
  Jet_X_Observables(type,xmin,xmax,nbins,mode,minn,maxn,lname,rname)
{
  m_name+="dR2_";
}

Primitive_Observable_Base * JetX_DeltaR_Distribution::Copy() const 
{
  JetX_DeltaR_Distribution * jdr =
    new JetX_DeltaR_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname,m_reflistname);
  return jdr;
}

double JetX_DeltaR_Distribution::Calc(const Vec4D &mom1,const Vec4D &mom2)
{
  double pt1  = mom1.PPerp();
  double pt2  = mom2.PPerp();
  double dphi = acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
  double deta = mom1.Eta()-mom2.Eta();
  return sqrt(sqr(deta) + sqr(dphi)); 
}

//----------------------------------------------------------------------

  class JetX_DeltaEta_Distribution : public Jet_X_Observables {
  protected:
  public:
    JetX_DeltaEta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
			      unsigned int mode, unsigned int minn, unsigned int maxn, 
			      const std::string & =std::string("Jets"),
			      const std::string & =std::string("FinalState"));

    Primitive_Observable_Base * Copy() const;
    double Calc(const Vec4D &m1,const Vec4D &m2);
  };

DEFINE_OBSERVABLE_GETTER(JetX_DeltaEta_Distribution,
			 JetX_DeltaEta_Distribution_Getter,"JetXDEta")

JetX_DeltaEta_Distribution::JetX_DeltaEta_Distribution(unsigned int type,double xmin,double xmax,int nbins,
							 unsigned int mode,unsigned int minn,unsigned int maxn, 
							 const std::string & lname,const std::string & rname) :
  Jet_X_Observables(type,xmin,xmax,nbins,mode,minn,maxn,lname,rname)
{
  m_name+="deta2_";
}

Primitive_Observable_Base * JetX_DeltaEta_Distribution::Copy() const 
{
  JetX_DeltaEta_Distribution * jde =
    new JetX_DeltaEta_Distribution(m_type,m_xmin,m_xmax,m_nbins,m_mode,m_minn,m_maxn,m_listname,m_reflistname);
  return jde;
}

double JetX_DeltaEta_Distribution::Calc(const Vec4D &mom1,const Vec4D &mom2)
{
  return dabs((mom1.Eta()-mom2.Eta()));
}
//----------------------------------------------------------------------

